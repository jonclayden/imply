#' Sparse images
#'
#' A sparse image stores only those spatial locations that hold data. Sparsity
#' is over *locations*, not over individual values: a location is either
#' present, in which case the whole vector of values held there is stored, or
#' absent, in which case every one of them is implicitly zero. That is what a
#' brain mask actually is, and it keeps the stored values contiguous.
#'
#' The mask is one bit per location and the values are packed in location
#' order, so finding a value costs a table lookup and a bit count rather than a
#' search. A coordinate list, the other common representation, would store
#' three or four indices alongside every value, and typically costs more memory
#' than it saves at the densities medical images actually have.
#'
#' A location holding `NA` is kept, since `NA` is not zero. Packing an image
#' therefore loses nothing.
#'
#' @param x An image, array, or for `asDense()` a sparse image.
#' @param ... Further arguments to `denseImage()` or `sparseImage()`.
#' @param mask A raw vector of one bit per location, or a logical vector.
#' @param values Packed values, in location order.
#' @param dim,spatial,pixdim,xform,spaceUnit,timeUnit Image geometry, as for
#'   [denseImage()].
#' @param template An image to take unspecified geometry from.
#' @name sparseImage
NULL

#' @rdname sparseImage
#' @export
sparseImage <- S7::new_class("sparseImage",
    properties = list(
        mask = S7::class_raw,
        values = S7::class_atomic,
        dims = S7::class_integer,
        spatial = S7::class_integer,
        pixdim = S7::class_double,
        xform = S7::class_double,
        spaceUnit = S7::class_character,
        timeUnit = S7::class_character
    ),
    validator = function (self) {
        nDims <- length(self@dims)

        if (length(self@spatial) != 1L || is.na(self@spatial))
            return("@spatial must be a single value")
        if (self@spatial < 0L || self@spatial > nDims)
            return(paste0("@spatial must be between 0 and ", nDims))
        if (anyNA(self@dims) || any(self@dims < 0L))
            return("@dims must not be missing or negative")

        if (!typeof(self@values) %in% c("logical", "integer", "double", "complex"))
            return("@values must be logical, integer, double or complex")

        locations <- prod(self@dims[seq_len(self@spatial)])
        if (length(self@mask) != ceiling(locations / 64) * 8)
            return("@mask is not the right length for the stated dimensions")

        elements <- prod(self@dims[-seq_len(self@spatial)])
        if (length(self@values) != maskCount(self@mask, locations) * elements)
            return("@values does not hold one entry per present location")

        if (length(self@pixdim) != self@spatial)
            return(paste0("@pixdim must have one element per spatial dimension (", self@spatial, ")"))
        if (!identical(dim(self@xform), c(4L, 4L)))
            return("@xform must be a 4x4 matrix")

        NULL
    },
    constructor = function (mask, values, dim, spatial = NULL, pixdim = NULL, xform = NULL,
                            spaceUnit = NULL, timeUnit = NULL, template = NULL)
    {
        dim <- as.integer(dim)
        nDims <- length(dim)
        spatial <- as.integer(spatial %||% attr(template, "spatial") %||% min(3L, nDims))

        if (is.logical(mask))
            mask <- maskFromLogical(mask)

        pixdim <- as.double(pixdim %||% attr(template, "pixdim") %||% rep(1, max(spatial, 0L)))
        xform <- xform %||% attr(template, "xform") %||% defaultXform(pixdim)
        xform <- as.matrix(xform)
        storage.mode(xform) <- "double"
        dimnames(xform) <- NULL

        S7::new_object(S7::S7_object(),
            mask = mask,
            values = values,
            dims = dim,
            spatial = spatial,
            pixdim = pixdim,
            xform = xform,
            spaceUnit = as.character(spaceUnit %||% attr(template, "spaceUnit") %||% "unknown"),
            timeUnit = as.character(timeUnit %||% attr(template, "timeUnit") %||% "unknown"))
    })

S7::S4_register(sparseImage)

#' @rdname sparseImage
#' @export
isSparseImage <- function (x) S7::S7_inherits(x, sparseImage)

#' @rdname sparseImage
#' @export
asSparse <- function (x, ...)
{
    if (isSparseImage(x))
        return(x)

    image <- asDenseImage(x, ...)
    packed <- denseToSparse(as.array(image), image@spatial)

    sparseImage(mask = packed$mask, values = packed$values, dim = dim(image),
                spatial = image@spatial, pixdim = image@pixdim, xform = image@xform,
                spaceUnit = image@spaceUnit, timeUnit = image@timeUnit)
}

#' @rdname sparseImage
#' @export
asDense <- function (x, ...)
{
    if (!isSparseImage(x))
        return(asDenseImage(x, ...))

    denseImage(sparseToDense(x@mask, x@values, x@dims, x@spatial),
               spatial = x@spatial, pixdim = x@pixdim, xform = x@xform,
               spaceUnit = x@spaceUnit, timeUnit = x@timeUnit)
}

#' @rdname sparseImage
#' @export
sparseness <- function (x)
{
    if (!isSparseImage(x))
        return(sum(as.array(x) == 0, na.rm = TRUE) / length(x))
    1 - maskCount(x@mask, locationCount(x)) / locationCount(x)
}

#' @rdname sparseImage
#' @export
mask <- function (x)
{
    if (!isSparseImage(x))
        stop("Only a sparse image carries a mask")
    array(maskToLogical(x@mask, locationCount(x)), x@dims[seq_len(x@spatial)])
}

locationCount <- function (x) prod(x@dims[seq_len(x@spatial)])
elementCount <- function (x) prod(x@dims[-seq_len(x@spatial)])

S7::method(dim, sparseImage) <- function (x) x@dims

S7::method(as.array, sparseImage) <- function (x, ...)
    sparseToDense(x@mask, x@values, x@dims, x@spatial)

S7::method(length, sparseImage) <- function (x) prod(x@dims)

S7::method(print, sparseImage) <- function (x, ...)
{
    present <- maskCount(x@mask, locationCount(x))

    cat(sprintf("Sparse image: %s (%s)\n", paste(x@dims, collapse = " x "), typeof(x@values)))
    if (x@spatial > 0L)
    {
        cat(sprintf("  Spatial dimensions : %s\n", paste(x@dims[seq_len(x@spatial)], collapse = " x ")))
        cat(sprintf("  Voxel dimensions   : %s %s\n",
                    paste(signif(x@pixdim, 4), collapse = " x "), x@spaceUnit))
        cat(sprintf("  Orientation        : %s\n", orientation(x)))
    }
    if (x@spatial < length(x@dims))
        cat(sprintf("  Values per location: %d\n", elementCount(x)))
    cat(sprintf("  Locations stored   : %s of %s (%.1f%% sparse)\n",
                format(present, big.mark = ","), format(locationCount(x), big.mark = ","),
                100 * sparseness(x)))

    invisible(x)
}

## Indexing goes through the mask rather than materialising the image, which
## is the whole point of the representation. As for dense images the result is
## a plain array: an arbitrary index has no well-defined geometry
S7::method(`[`, sparseImage) <- function (x, ..., drop = TRUE)
{
    ## x[] supplies one argument which is nonetheless missing, so the empty
    ## symbol has to be looked for rather than the arguments merely counted
    indices <- as.list(substitute(list(...)))[-1L]
    supplied <- length(indices)
    absent <- vapply(indices, identical, NA, quote(expr = ))

    if (supplied == 0L || all(absent))
        return(as.array(x))

    if (supplied == 1L)
    {
        i <- ..1
        if (is.matrix(i) && ncol(i) == length(x@dims))
            i <- flattenIndices(array(0L, x@dims), i)
        else if (is.logical(i))
            i <- which(i)
        return(sparseElements(x@mask, x@values, x@dims, x@spatial, as.double(i)))
    }

    ## Full n-dimensional indexing is rare enough on a sparse image that
    ## materialising is the simpler and safer path
    call <- sys.call()
    call[[1L]] <- quote(`[`)
    call[[2L]] <- as.array(x)
    eval(call, parent.frame())
}
