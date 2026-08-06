#' Dense images
#'
#' A dense image is an array carrying image geometry alongside it. It is an S7
#' class whose parent is the `array` S3 class, which means every storage mode
#' works (double, integer, logical and complex all occur in practice, masks
#' being logical), the data are still a plain R array so no copy is needed to
#' pass them to compiled code, and the geometry is validated whenever it is
#' set rather than only on construction.
#'
#' Properties are stored as ordinary attributes, so compiled code reads them
#' without needing to know anything about S7.
#'
#' @param .data An array, or any atomic vector, which is treated as
#'   one-dimensional.
#' @param pixdim Voxel dimensions, one per spatial dimension.
#' @param xform A 4x4 affine transform mapping voxel to world coordinates.
#' @param spatial The number of leading dimensions that index location rather
#'   than the value held at each location. Defaults to three, or the
#'   dimensionality if that is smaller.
#' @param spaceUnit,timeUnit Units of measurement.
#' @param template An image to take unspecified geometry from.
#' @param x An image, or an object to coerce to one.
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
        pixdim = S7::class_double,
        xform = S7::class_double,
        spaceUnit = S7::class_character,
        timeUnit = S7::class_character
    ),
    validator = function (self) {
        nDims <- length(dim(self))

        if (length(self@spatial) != 1L || is.na(self@spatial))
            return("@spatial must be a single value")
        if (self@spatial < 0L || self@spatial > nDims)
            return(paste0("@spatial must be between 0 and ", nDims))

        if (length(self@pixdim) != self@spatial)
            return(paste0("@pixdim must have one element per spatial dimension (", self@spatial, ")"))
        if (anyNA(self@pixdim))
            return("@pixdim must not be missing")

        if (!identical(dim(self@xform), c(4L, 4L)))
            return("@xform must be a 4x4 matrix")
        if (anyNA(self@xform))
            return("@xform must not contain missing values")
        if (!isTRUE(all.equal(self@xform[4, ], c(0, 0, 0, 1))))
            return("@xform must be affine, with a final row of (0, 0, 0, 1)")

        if (length(self@spaceUnit) != 1L || length(self@timeUnit) != 1L)
            return("@spaceUnit and @timeUnit must each be a single value")

        NULL
    },
    constructor = function (.data, pixdim = NULL, xform = NULL, spatial = NULL,
                            spaceUnit = NULL, timeUnit = NULL, template = NULL)
    {
        if (!is.atomic(.data))
            stop("Image data must be an atomic array")
        if (!typeof(.data) %in% c("logical", "integer", "double", "complex"))
            stop("Image data must be logical, integer, double or complex, not ", typeof(.data))
        if (is.null(dim(.data)))
            dim(.data) <- length(.data)

        nDims <- length(dim(.data))

        ## Explicit arguments win, then the template, then defaults
        spatial <- as.integer(spatial %||% attr(template, "spatial") %||% min(3L, nDims))
        pixdim <- as.double(pixdim %||% attr(template, "pixdim") %||% rep(1, max(spatial, 0L)))
        xform <- xform %||% attr(template, "xform") %||% defaultXform(pixdim)

        xform <- as.matrix(xform)
        storage.mode(xform) <- "double"
        dimnames(xform) <- NULL

        S7::new_object(.data,
            spatial = spatial,
            pixdim = pixdim,
            xform = xform,
            spaceUnit = as.character(spaceUnit %||% attr(template, "spaceUnit") %||% "unknown"),
            timeUnit = as.character(timeUnit %||% attr(template, "timeUnit") %||% "unknown"))
    })

#' @rdname denseImage
#' @export
#' @details
#' Note that S7 qualifies a class name with its package, so the class
#' attribute is `"imply::denseImage"` and `inherits(x, "denseImage")` is
#' `FALSE`. Use `isDenseImage()` rather than testing the class directly.
isDenseImage <- function (x) S7::S7_inherits(x, denseImage)

#' @rdname denseImage
#' @export
asDenseImage <- function (x, ...)
{
    if (isDenseImage(x))
        x
    else
        denseImage(x, ...)
}

## Subsetting. S7 objects are not subsettable by default, so without these an
## image could not be indexed at all.
##
## The result is a plain array rather than an image: an arbitrary index has no
## well-defined geometry, and this matches both base R's behaviour for a
## classed array and the convention established by tractor.base. Use crop()
## when the geometry should be carried through.

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
        cat(sprintf("  Voxel dimensions   : %s %s\n",
                    paste(signif(x@pixdim, 4), collapse = " x "), x@spaceUnit))
        cat(sprintf("  Orientation        : %s\n", orientation(x)))
    }
    if (nSpatial < length(dims))
        cat(sprintf("  Values per location: %d\n", prod(dims[-seq_len(nSpatial)])))

    invisible(x)
}

S7::method(as.array, denseImage) <- function (x, ...)
{
    for (name in c("spatial", "pixdim", "xform", "spaceUnit", "timeUnit", "S7_class"))
        attr(x, name) <- NULL
    class(x) <- NULL
    x
}
