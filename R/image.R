#' Dense images
#'
#' A dense image is an array carrying image geometry alongside it. It is an S7
#' class whose parent is the `array` S3 class, which means every storage mode
#' works (double, integer, logical and complex all occur in practice, masks
#' being logical). The data are still a plain R array, so no copy is needed to
#' pass them to compiled code, and the geometry is validated whenever it is
#' set rather than only on construction.
#'
#' Properties are stored as ordinary attributes, so compiled code reads them
#' without needing to know anything about S7.
#' 
#' @note S7 qualifies a class name with its package, so the class
#' attribute is `"imply::denseImage"` and `inherits(x, "denseImage")` is
#' `FALSE`. Use `isDenseImage()` rather than testing the class directly.
#'
#' @param .data An array, or any atomic vector, which is treated as
#'   one-dimensional.
#' @param voxelSize Voxel size, one per spatial dimension.
#' @param worldTransform A 4x4 affine transform mapping voxel to world
#'   coordinates. Decomposed into rotation/translation and voxel size on
#'   assignment; see [geometry].
#' @param spatial The number of leading dimensions that index location rather
#'   than the value held at each location. Defaults to three, or the
#'   dimensionality if that is smaller.
#' @param spaceUnit, timeUnit Units of measurement.
#' @param template An image to take unspecified geometry from.
#' @param x An image.
#' @param ... Further arguments to `denseImage()`.
#' @name denseImage
NULL

## S7 needs to know how to build the parent, which for an array just means
## ensuring there is a dim attribute
arrayClass <- S7::new_S3_class("array", constructor = function (.data = array(numeric())) {
    if (is.null(dim(.data)))
        dim(.data) <- length(.data)
    .data
})

#' @rdname denseImage
#' @export
denseImage <- S7::new_class("denseImage",
    parent = arrayClass,
    properties = list(
        spatial = S7::class_integer,
        voxelSize = S7::class_double,
        orientation = S7::class_double,
        spaceUnit = S7::class_character,
        timeUnit = S7::class_character
    ),
    validator = function (self) {
        nDims <- length(dim(self))

        if (length(self@spatial) != 1L || is.na(self@spatial))
            return("@spatial must be a single value")
        if (self@spatial < 0L || self@spatial > nDims)
            return(paste0("@spatial must be between 0 and ", nDims))

        if (length(self@voxelSize) != self@spatial)
            return(paste0("@voxelSize must have one element per spatial dimension (", self@spatial, ")"))
        if (anyNA(self@voxelSize))
            return("@voxelSize must not be missing")
        if (any(self@voxelSize <= 0))
            return("@voxelSize must be strictly positive")

        if (!identical(dim(self@orientation), c(4L, 4L)))
            return("@orientation must be a 4x4 matrix")
        if (anyNA(self@orientation))
            return("@orientation must not contain missing values")
        if (!isTRUE(all.equal(self@orientation[4, ], c(0, 0, 0, 1))))
            return("@orientation must be affine, with a final row of (0, 0, 0, 1)")
        block <- self@orientation[1:3, 1:3, drop = FALSE]
        if (max(abs(crossprod(block) - diag(3))) > orthogonalityTolerance)
            return("@orientation must be rigid: a rotation or reflection, with no scale or shear")

        if (length(self@spaceUnit) != 1L || length(self@timeUnit) != 1L)
            return("@spaceUnit and @timeUnit must each be a single value")

        NULL
    },
    constructor = function (.data, voxelSize = NULL, worldTransform = NULL, spatial = NULL,
                            spaceUnit = NULL, timeUnit = NULL, template = NULL)
    {
        if (!is.atomic(.data))
            stop("Image data must be an atomic array")
        if (!typeof(.data) %in% c("logical", "integer", "double", "complex"))
            stop("Image data must be logical, integer, double or complex, not ", typeof(.data))
        if (is.null(dim(.data)))
            dim(.data) <- length(.data)

        nDims <- length(dim(.data))

        ## Explicit arguments win, then a value implied by another explicit
        ## argument (a worldTransform implies both orientation and voxel
        ## size), then the template, then defaults
        spatial <- as.integer(spatial %||% attr(template, "spatial") %||% min(3L, nDims))

        decomposed <- if (is.null(worldTransform)) NULL
                      else decomposeTransform(validateXform(worldTransform), spatial)
        orientation <- decomposed$orientation %||% attr(template, "orientation") %||% diag(4)
        voxelSize <- as.double(voxelSize %||% decomposed$voxelSize %||%
                               attr(template, "voxelSize") %||% rep(1, max(spatial, 0L)))

        S7::new_object(.data,
            spatial = spatial,
            voxelSize = voxelSize,
            orientation = orientation,
            spaceUnit = as.character(spaceUnit %||% attr(template, "spaceUnit") %||% "unknown"),
            timeUnit = as.character(timeUnit %||% attr(template, "timeUnit") %||% "unknown"))
    })

#' @rdname denseImage
#' @export
isDenseImage <- function (x) S7::S7_inherits(x, denseImage)

#' @rdname denseImage
#' @export
asDense <- function (x, ...)
{
    ## Unpacks whichever of the compact representations it is given, so a
    ## caller that just wants ordinary values need not ask which one it has.
    if (isDenseImage(x))
        return (x)
    else if (isPackedImage(x))
        return(denseImage(as.array(x), spatial = x@spatial, voxelSize = x@voxelSize,
                          worldTransform = worldTransform(x),
                          spaceUnit = x@spaceUnit, timeUnit = x@timeUnit))
    else if (isSparseImage(x))
        return(denseImage(sparseToDense(x@mask, x@values, x@dims, x@spatial),
                          spatial = x@spatial, voxelSize = x@voxelSize,
                          worldTransform = worldTransform(x), spaceUnit = x@spaceUnit,
                          timeUnit = x@timeUnit))
    else
        return(denseImage(x, ...))
}

## S7 objects are not subsettable by default, so without these methods an image
## could not be indexed at all
##
## The result is a plain array rather than an image: an arbitrary index has no
## well-defined geometry, and this matches both base R's behaviour for a
## classed array and the convention established by tractor.base. Use crop()
## when the geometry should be carried through
##
## The call is rebuilt rather than forwarded, because missing index arguments
## (as in x[,,1]) cannot be passed through `...` faithfully
S7::method(`[`, denseImage) <- function (x, ..., drop = TRUE)
{
    call <- sys.call()
    call[[1L]] <- quote(`[`)
    call[[2L]] <- as.array(x)
    eval(call, parent.frame())
}

S7::method(`[<-`, denseImage) <- function (x, ..., value)
{
    call <- sys.call()
    call[[1L]] <- quote(`[<-`)
    call[[2L]] <- as.array(x)
    result <- eval(call, parent.frame())

    ## Replacement preserves shape, so the geometry still applies
    denseImage(result, template = x)
}

S7::method(print, denseImage) <- function (x, ...)
{
    dims <- dim(x)
    nSpatial <- x@spatial

    cat(sprintf("Dense image: %s (%s)\n", paste(dims, collapse = " x "), typeof(x)))
    if (nSpatial > 0L)
    {
        cat(sprintf("  Spatial dimensions : %s\n", paste(dims[seq_len(nSpatial)], collapse = " x ")))
        cat(sprintf("  Voxel size         : %s %s\n",
                    paste(signif(x@voxelSize, 4), collapse = " x "), x@spaceUnit))
    }
    if (nSpatial < length(dims))
        cat(sprintf("  Values per location: %d\n", prod(dims[-seq_len(nSpatial)])))

    invisible(x)
}

S7::method(as.array, denseImage) <- function (x, ...)
{
    for (name in c("spatial", "voxelSize", "orientation", "spaceUnit", "timeUnit", "S7_class"))
        attr(x, name) <- NULL
    class(x) <- NULL
    x
}
