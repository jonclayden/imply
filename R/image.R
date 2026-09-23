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
#'   than the value held at each location. Defaults to the number the geometry
#'   has, if one is given, or otherwise to three, or the dimensionality if that
#'   is smaller.
#' @param unit The unit of measurement of voxel size and world coordinates.
#' @param geometry An [imageGeometry][geometry], or an image to take one from, giving
#'   whatever is not specified explicitly. Its grid must match the leading
#'   dimensions of the data.
#' @param x An image.
#' @param ... Further arguments to `denseImage()`.
#' @return An object of S7 class `denseImage` representing a dense image, with
#'   a `geometry` property holding its [imageGeometry][geometry].
#' @name denseImage
NULL

## S7 needs to know how to build the parent, which for an array just means
## ensuring there is a dim attribute
arrayClass <- S7::new_S3_class("array", constructor = function (.data = array(numeric())) {
    if (is.null(dim(.data)))
        dim(.data) <- length(.data)
    .data
})

## The geometry must describe the leading dimensions of the data, which is the
## one thing about it that the geometry cannot check for itself
geometryMismatch <- function (geometry, dims)
{
    spatial <- length(geometry@dims)
    if (spatial > length(dims))
        return(paste0("@geometry has ", spatial, " spatial dimensions, but the image has only ", length(dims)))
    if (!identical(geometry@dims, as.integer(dims[seq_len(spatial)])))
        return(paste0("@geometry is for a grid of ", formatDims(geometry@dims),
                      ", which does not match the image's leading dimensions"))
    NULL
}

#' @rdname denseImage
#' @export
denseImage <- S7::new_class("denseImage",
    parent = arrayClass,
    properties = list(
        geometry = imageGeometry
    ),
    validator = function (self) geometryMismatch(self@geometry, dim(self)),
    constructor = function (.data, voxelSize = NULL, worldTransform = NULL, spatial = NULL,
                            unit = NULL, geometry = NULL)
    {
        if (!is.atomic(.data))
            stop("Image data must be an atomic array")
        if (!typeof(.data) %in% c("logical", "integer", "double", "complex"))
            stop("Image data must be logical, integer, double or complex, not ", typeof(.data))
        if (is.null(dim(.data)))
            dim(.data) <- length(.data)

        S7::new_object(.data,
            geometry = resolveGeometry(dim(.data), spatial, voxelSize, worldTransform, unit, geometry))
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
        return(denseImage(as.array(x), geometry = x@geometry))
    else if (isSparseImage(x))
        return(denseImage(sparseToDense(x@mask, x@values, x@dims, spatial(x)), geometry = x@geometry))
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

## Registered as an ordinary call rather than via `[<-`(x) <- value sugar,
## because that sugar assigns its result back into a package-namespace
## binding literally named `[<-`. R CMD check's "checking replacement
## functions" step then tries to inspect that binding and crashes on some
## older R versions (get(f, envir = code_env) : invalid first argument),
## rather than just noting it. method<- is an ordinary function, so calling
## it directly avoids creating the binding at all
S7::`method<-`(`[<-`, denseImage, function (x, ..., value)
{
    call <- sys.call()
    call[[1L]] <- quote(`[<-`)
    call[[2L]] <- as.array(x)
    result <- eval(call, parent.frame())

    ## Replacement preserves shape, so the geometry still applies
    denseImage(result, geometry = x@geometry)
})

S7::method(print, denseImage) <- function (x, ...)
{
    dims <- dim(x)
    nSpatial <- spatial(x)

    cat(sprintf("Dense image: %s (%s)\n", paste(dims, collapse = " x "), typeof(x)))
    printGeometry(x@geometry)
    if (nSpatial < length(dims))
        cat(sprintf("  Values per location: %d\n", prod(dims[-seq_len(nSpatial)])))

    invisible(x)
}

S7::method(as.array, denseImage) <- function (x, ...)
{
    for (name in c("geometry", "S7_class"))
        attr(x, name) <- NULL
    class(x) <- NULL
    x
}
