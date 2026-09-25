#' Sparse images
#'
#' A sparse image stores only those spatial locations that hold data. Sparsity
#' is over *locations*, not over individual values: a location is either
#' present, in which case the whole vector of values held there is stored, or
#' absent, in which case every one of them is implicitly zero. That generally
#' matches the intent, and it keeps the stored values contiguous.
#'
#' The mask is one bit per location and the values are packed in location
#' order, so finding a value costs a table lookup and a bit count rather than a
#' search. A coordinate list, the other common representation, would store
#' three or four indices alongside every value, and may cost more memory than
#' saves at typical densities.
#'
#' A location holding `NA` is kept, since `NA` is not zero. Packing an image
#' therefore loses nothing.
#'
#' A sparse image can be indexed and assigned to like an array. Assignment
#' keeps the mask accurate: a location given a non-zero value is added, and one
#' left holding only zeros is dropped. That takes a pass over the stored
#' values each time, which suits occasional edits; to set many values, gather
#' them and assign once, rather than one at a time in a loop. As for an array,
#' values are promoted to a wider replacement type.
#'
#' @param x An image or array.
#' @param ... Further arguments to `sparseImage()`.
#' @param mask A raw vector of one bit per location, or a logical vector.
#' @param values Packed values, in location order.
#' @param dim The full dimensions of the image.
#' @param spatial,voxelSize,worldTransform,unit,geometry Image geometry, as
#'   for [denseImage()].
#' @param layout How the image's axes map onto the stored order of locations,
#'   as described for [storageLayout()]. For a sparse image it may only
#'   reorder and reverse the spatial axes among themselves; the values held at
#'   each location are always stored in order. The default is the identity.
#' @return An object of S7 class `sparseImage` representing a sparse image,
#'   with properties corresponding to the arguments listed above.
#' @name sparseImage
NULL

#' @rdname sparseImage
#' @export
sparseImage <- S7::new_class("sparseImage",
    properties = list(
        mask = S7::class_raw,
        values = S7::class_atomic,
        dims = S7::class_integer,
        layout = S7::class_integer,
        geometry = imageGeometry
    ),
    validator = function (self) {
        if (anyNA(self@dims) || any(self@dims < 0L))
            return("@dims must not be missing or negative")
        mismatch <- geometryMismatch(self@geometry, self@dims) %||% checkLayout(self@layout, length(self@dims))
        if (!is.null(mismatch))
            return(mismatch)

        ## The accessor splits a stored index into a location and an element,
        ## which only works while the spatial axes stay in front
        nSpatial <- spatial(self)
        if (any(abs(self@layout[seq_len(nSpatial)]) > nSpatial) ||
            !isIdentityLayout(self@layout[-seq_len(nSpatial)] - nSpatial))
            return("@layout of a sparse image may only reorder its spatial axes among themselves")

        if (!typeof(self@values) %in% c("logical", "integer", "double", "complex"))
            return("@values must be logical, integer, double or complex")

        locations <- locationCount(self)
        if (length(self@mask) != ceiling(locations / 64) * 8)
            return("@mask is not the right length for the stated dimensions")

        if (length(self@values) != maskCount(self@mask, locations) * elementCount(self))
            return("@values does not hold one entry per present location")

        NULL
    },
    constructor = function (mask, values, dim, spatial = NULL, voxelSize = NULL, worldTransform = NULL,
                            unit = NULL, geometry = NULL, layout = NULL)
    {
        dim <- as.integer(dim)
        nDims <- length(dim)
        geometry <- resolveGeometry(dim, spatial, voxelSize, worldTransform, unit, geometry)
        spatial <- length(geometry@dims)

        if (is.logical(mask))
            mask <- maskFromLogical(mask)

        ## Values are stored already shaped, one column per stored location,
        ## so that maskedMatrix() can hand them back without copying. R would
        ## otherwise duplicate them on the way out, since they are shared with
        ## this object
        elements <- if (spatial < nDims) prod(dim[-seq_len(spatial)]) else 1L
        if (length(values) > 0L || elements > 0L)
            dim(values) <- c(elements, length(values) %/% max(elements, 1L))

        S7::new_object(S7::S7_object(),
            mask = mask,
            values = values,
            dims = dim,
            layout = as.integer(layout %||% seq_along(dim)),
            geometry = geometry)
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

    image <- asDense(x, ...)
    packed <- denseToSparse(as.array(image), spatial(image))
    sparseImage(mask = packed$mask, values = packed$values, dim = dim(image), geometry = image@geometry)
}

#' `maskedMatrix()` returns the stored values, with one column per stored
#' location, which is the data matrix most voxelwise analysis wants: nothing
#' outside the mask is present at all, and the values at one location are
#' contiguous. The columns are in the order of `which(mask(x))`. Packing is
#' voxel-major, so this is the shape the values are already held in, and
#' returning them costs nothing unless the image has been reoriented, in which
#' case the columns are reordered to match, which copies.
#'
#' `storedValues()` returns the values exactly as stored, never copying, along
#' with the location of each column as a linear index into the image's own
#' spatial grid. The two differ only for a reoriented image.
#'
#' @rdname sparseImage
#' @export
maskedMatrix <- function (x)
{
    if (!isSparseImage(x))
        stop("Only a sparse image has packed values")

    if (!isIdentityLayout(x@layout))
    {
        stored <- storedValues(x)
        return(stored$values[, order(stored$locations), drop = FALSE])
    }

    storedMatrix(x)
}

#' @rdname sparseImage
#' @export
storedValues <- function (x)
{
    if (!isSparseImage(x))
        stop("Only a sparse image has packed values")

    present <- which(maskToLogical(x@mask, locationCount(x)))
    locations <- if (isIdentityLayout(x@layout)) present else {
        ## Where each view location is stored, inverted to find where each
        ## stored location sits in the view
        nSpatial <- spatial(x)
        viewToStored <- viewIndices(seq_len(locationCount(x)), x@dims[seq_len(nSpatial)], x@layout[seq_len(nSpatial)])
        storedToView <- integer(length(viewToStored))
        storedToView[viewToStored] <- seq_along(viewToStored)
        storedToView[present]
    }

    list(values = storedMatrix(x), locations = locations)
}

## The values as held, shaped with one column per stored location
storedMatrix <- function (x)
{
    values <- x@values
    elements <- elementCount(x)
    stored <- maskCount(x@mask, locationCount(x))

    ## Already shaped by the constructor; reshaped only for an image built
    ## some other way, which is the one case that has to copy
    if (identical(dim(values), c(as.integer(elements), as.integer(stored))))
        return(values)

    dim(values) <- c(elements, stored)
    values
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
    present <- maskToLogical(x@mask, locationCount(x))
    nSpatial <- spatial(x)
    if (isIdentityLayout(x@layout))
        array(present, x@geometry@dims)
    else
        viewGather(present, x@geometry@dims, x@layout[seq_len(nSpatial)])
}

#' Masks over spatial locations
#'
#' A mask selects spatial locations, for example to restrict [voxelApply()] to
#' a region of interest. It may be given in several forms, and this function
#' reduces any of them to one logical value per location.
#'
#' @param mask A logical array, a numeric array in which non-zero values are
#'   selected, a sparse image whose own mask is used, or any other image,
#'   which is treated as a numeric array.
#' @param spatialDims The spatial dimensions the mask must cover, or an image
#'   or [imageGeometry][geometry] to take them from.
#' @return A logical vector with one element per spatial location.
#' @export
asMaskVector <- function (mask, spatialDims)
{
    if (isImage(spatialDims) || isImageGeometry(spatialDims))
        spatialDims <- geometry(spatialDims)@dims

    if (isSparseImage(mask))
        mask <- mask(mask)
    else if (isImage(mask))
        mask <- as.array(asDense(mask))

    if (is.logical(mask))
        selected <- as.vector(mask)
    else if (is.numeric(mask))
        selected <- as.vector(mask) != 0
    else
        stop("Mask must be a logical array, a numeric array, or an image")

    if (anyNA(selected))
        stop("Mask must not contain missing values")
    if (length(selected) != prod(spatialDims))
        stop("Mask must have one element per spatial location (", prod(spatialDims), ")")

    selected
}

locationCount <- function (x) prod(x@geometry@dims)
elementCount <- function (x) prod(x@dims[-seq_len(spatial(x))])

S7::method(dim, sparseImage) <- function (x) x@dims

S7::method(as.array, sparseImage) <- function (x, ...)
    sparseToDense(x@mask, x@values, x@dims, spatial(x), layoutArg(x))

S7::method(length, sparseImage) <- function (x) prod(x@dims)

S7::method(print, sparseImage) <- function (x, ...)
{
    present <- maskCount(x@mask, locationCount(x))

    cat(sprintf("Sparse image: %s (%s)\n", paste(x@dims, collapse = " x "), typeof(x@values)))
    printGeometry(x@geometry)
    if (spatial(x) < length(x@dims))
        cat(sprintf("  Values per location: %d\n", elementCount(x)))
    cat(sprintf("  Locations stored   : %s of %s (%.1f%% sparse)\n",
                format(present, big.mark = ","), format(locationCount(x), big.mark = ","),
                100 * sparseness(x)))
    printLayout(x@layout)

    invisible(x)
}

