#' Apply a function over an image
#'
#' `imapply()` is a memory-efficient analogue of [base::apply()]. It gives the
#' same answer, but where `apply()` permutes the whole array into a fresh copy
#' before looping, `imapply()` gathers each sub-array directly through the
#' stride vector. Peak memory is therefore the input plus the result, rather
#' than twice the input plus the result, which matters once an image is larger
#' than a comfortable fraction of memory.
#'
#' The remaining functions differ only in how much of the *space* they hand to
#' `fun` at a time: a single location for `voxelApply()`, a one-dimensional
#' line for `lineApply()`, and a two-dimensional slice for `sliceApply()`.
#'
#' In every case the values held at each location travel with the unit. So for
#' an image with a time series at each location, `voxelApply()` passes one
#' series, `lineApply()` passes a line's worth of series, and `sliceApply()`
#' passes a slice's worth:
#'
#' \preformatted{
#'   dim(x) = 4 x 5 x 6 x 10, spatial = 3
#'
#'   voxelApply(x, f)           f sees  10 values
#'   lineApply(x, f, axis = 1)  f sees  4 x 10
#'   sliceApply(x, f, axis = 3) f sees  4 x 5 x 10
#' }
#'
#' `axis` always names a spatial dimension, but means what the name of each
#' function implies: a line *runs along* its axis, whereas a slice is *cut
#' across* its axis. Use `imapply()` directly to iterate over a non-spatial
#' dimension, such as applying a function to each volume in a time series.
#'
#' @param x An image or a plain array.
#' @param margin The dimensions to retain, as for [base::apply()].
#' @param fun A function to apply.
#' @param ... Further arguments to `fun`.
#' @param simplify Whether to simplify the result to an array where possible.
#' @param threads Number of threads to use, or `NULL` to consult
#'   `getOption("imply.threads")`. See [parallelism].
#' @param progress `FALSE` for none, `TRUE` for a text progress bar showing
#'   the percentage complete and the rate in voxels per second, or a function
#'   of `(done, total)`.
#' @param axis For `lineApply()`, the axis lines run along; for `sliceApply()`,
#'   the axis slices cut across.
#' @param mask For `voxelApply()`, a logical array over the spatial dimensions,
#'   a sparse image whose mask is to be used, or `NULL` for none. Locations
#'   outside it are not visited at all.
#' @param fill The value given to locations outside `mask`.
#' @return For `imapply()`, as [base::apply()]. For `voxelApply()`, an image
#'   when the function returns a single value per location, otherwise an array.
#' @name imapply
NULL

#' @rdname imapply
#' @export
imapply <- function (x, margin, fun, ..., simplify = TRUE, threads = NULL, progress = FALSE)
{
    fun <- match.fun(fun)

    dims <- dim(x)
    if (is.null(dims))
        dims <- length(x)
    margin <- as.integer(margin)

    if (length(margin) == 0L)
        stop("Margin must name at least one dimension")
    if (anyNA(margin) || any(margin < 1L) || any(margin > length(dims)))
        stop("Margin is out of range for an array with ", length(dims), " dimensions")
    if (anyDuplicated(margin))
        stop("Margin contains a repeated dimension")

    ## Wrapping here means the compiled loop only ever calls a function of one
    ## argument, so the call can be built once and reused
    wrapped <- function (v) fun(v, ...)

    marginDims <- dims[margin]
    callDims <- dims[-margin]
    nCalls <- prod(marginDims)

    dimNames <- dimnames(x)
    marginNames <- if (!is.null(dimNames)) dimNames[margin] else NULL
    callNames <- if (!is.null(dimNames)) dimNames[-margin] else NULL
    if (isAllNull(marginNames))
        marginNames <- NULL
    if (isAllNull(callNames))
        callNames <- NULL

    ## Tested before routing, because base::apply() calls the function once on
    ## a dummy sub-array to establish the type of an empty result, and the
    ## compiled kernels have nothing to say about that
    if (nCalls == 0L)
        return(applyToEmpty(x, wrapped, callDims, callNames, marginDims, marginNames))

    ## A recognised reduction is answered by the compiled kernel instead, which
    ## avoids an interpreter call per sub-array and can use worker threads. The
    ## check is conservative, so this never changes the answer
    extras <- list(...)
    routed <- routableReduction(fun, x, extras, simplify)
    if (!is.null(routed))
        return(imreduce(x, margin, routed, na.rm = isTRUE(extras$na.rm), threads = threads))

    reporter <- newProgress(progress, nCalls, progressUnit(x, dims, margin))
    if (!is.null(reporter))
        on.exit(reporter$close(), add = TRUE)

    out <- runOverMargin(unclassArray(x), margin, wrapped, callNames, simplify,
                         nCalls, resolveThreads(threads), reporter)
    shapeResult(out, marginDims, marginNames, margin, simplify)
}

isAllNull <- function (x) is.null(x) || all(vapply(x, is.null, logical(1)))

## The result-shaping rules are base::apply()'s, reproduced rather than
## reinvented so that the two agree exactly. In particular a run of calls that
## all return NULL simplifies to NULL, not to a list of NULLs
shapeResult <- function (out, marginDims, marginNames, margin, simplify)
{
    nCalls <- prod(marginDims)
    innerNames <- NULL

    if (out$isList)
    {
        results <- out$values
        first <- results[[1L]]

        asList <- !simplify || is.recursive(first)
        innerLength <- length(first)
        innerNames <- names(first)

        if (!asList)
            asList <- any(lengths(results) != innerLength)
        if (!asList && length(innerNames) > 0L)
        {
            same <- vapply(results, function (r) identical(names(r), innerNames), NA)
            if (!all(same))
                innerNames <- NULL
        }

        ## A list that cannot be simplified is still shaped over the margins,
        ## which for more than one margin means an array of mode list
        if (asList)
        {
            values <- results
            total <- nCalls
            innerNames <- NULL
        }
        else
        {
            values <- unlist(results, recursive = FALSE)
            total <- length(values)
        }
    }
    else
    {
        values <- out$values
        total <- length(values)
    }

    if (length(margin) == 1L && total == nCalls)
    {
        if (!is.null(marginNames))
            names(values) <- marginNames[[1L]]
        values
    }
    else if (total == nCalls)
        array(values, marginDims, marginNames)
    else if (total > 0L && total %% nCalls == 0L)
    {
        allNames <- if (is.null(innerNames) && is.null(marginNames)) NULL
                    else c(list(innerNames),
                           if (is.null(marginNames)) rep(list(NULL), length(marginDims)) else marginNames)
        array(values, c(total %/% nCalls, marginDims), allNames)
    }
    else
        values
}

## When a margin has extent zero there is nothing to iterate over, but the
## function is still called once on a dummy sub-array to establish the type of
## the empty result
applyToEmpty <- function (x, wrapped, callDims, callNames, marginDims, marginNames)
{
    dummy <- array(vector(typeof(x), 1L), dim = c(prod(callDims), 1L))
    value <- wrapped(if (length(callDims) < 2L) dummy[, 1L] else array(dummy[, 1L], callDims, callNames))

    if (is.null(value))
        value
    else if (length(marginDims) < 2L)
        value[1L][-1L]
    else
        array(value, marginDims, marginNames)
}

## Packed and sparse images are passed through untouched, since the compiled
## loop reads them in place. Only a dense image needs unclassing, which avoids
## any chance of a method being invoked during the loop
unclassArray <- function (x)
{
    if (isPackedImage(x) || isSparseImage(x))
        x
    else if (isDenseImage(x))
        as.array(x)
    else
        x
}

#' @rdname imapply
#' @export
voxelApply <- function (x, fun, ..., mask = NULL, fill = 0, simplify = TRUE,
                        threads = NULL, progress = FALSE)
{
    nSpatial <- spatial(x)
    dims <- dim(x)
    if (is.null(dims))
        dims <- length(x)

    if (nSpatial < 1L)
        stop("Image has no spatial dimensions to apply over")
    if (nSpatial == length(dims))
        stop("Image holds a single value at each location, so there is nothing to apply over")

    if (!is.null(mask))
        return(maskedVoxelApply(x, fun, ..., mask = mask, fill = fill, simplify = simplify,
                                threads = threads, progress = progress))

    result <- imapply(x, seq_len(nSpatial), fun, ..., simplify = simplify, threads = threads, progress = progress)

    ## A single value per location is itself an image, and inherits the
    ## geometry of the input, whichever way that input was stored
    if (isImage(x) && simplify && is.atomic(result) &&
        typeof(result) %in% c("logical", "integer", "double", "complex") &&
        length(result) == prod(dims[seq_len(nSpatial)]))
        result <- denseImage(array(result, dims[seq_len(nSpatial)]), template = x, spatial = nSpatial)

    result
}

#' @rdname imapply
#' @export
lineApply <- function (x, fun, ..., axis = 1L, simplify = TRUE, threads = NULL, progress = FALSE)
{
    nSpatial <- spatial(x)
    axis <- checkAxis(axis, nSpatial)

    if (nSpatial < 2L)
        stop("An image with one spatial dimension is a single line, so apply the function to it directly")

    ## Every spatial dimension but the axis is retained, so fun sees one line
    ## running along the axis, together with the values at each of its
    ## locations. Lines never overlap, which is what makes this decomposition
    ## safe to parallelise
    imapply(x, seq_len(nSpatial)[-axis], fun, ..., simplify = simplify, threads = threads, progress = progress)
}

checkAxis <- function (axis, nSpatial)
{
    axis <- as.integer(axis)
    if (length(axis) != 1L || is.na(axis) || axis < 1L || axis > nSpatial)
        stop("Axis must be a single spatial dimension of the image (1 to ", nSpatial,
             "); use imapply() to iterate over a non-spatial dimension")
    axis
}

#' @rdname imapply
#' @export
sliceApply <- function (x, fun, ..., axis = 3L, simplify = TRUE, threads = NULL, progress = FALSE)
{
    nSpatial <- spatial(x)

    if (nSpatial < 3L)
        stop("Slices need three spatial dimensions; use lineApply() for a two-dimensional image")
    axis <- checkAxis(axis, nSpatial)

    ## Only the axis is retained, so fun sees the plane cut across it, together
    ## with the values at each of its locations
    imapply(x, axis, fun, ..., simplify = simplify, threads = threads, progress = progress)
}


## Applying only where a mask holds
##
## Nothing is gained by walking a sparse image and gathering zeros for the
## locations it does not store: that still makes one call per location. What
## saves the work is doing the loop in the packed space, over a matrix with one
## column per selected location, and scattering the answers back afterwards
##
## When the image is already sparse and the mask is its own, that matrix is the
## stored values themselves and costs nothing. Otherwise the selected values
## are gathered once
maskedVoxelApply <- function (x, fun, ..., mask, fill, simplify, threads, progress)
{
    dims <- dim(x)
    if (is.null(dims))
        dims <- length(x)
    nSpatial <- spatial(x)
    spatialDims <- dims[seq_len(nSpatial)]
    nLocations <- prod(spatialDims)

    selected <- asMaskVector(mask, spatialDims)
    index <- which(selected)

    if (length(index) == 0L)
        stop("Mask selects no locations")

    ## The packed matrix is not an image, so the unit imapply() would infer for
    ## the rate is not the right one. Here one call is exactly one location,
    ## whatever shape the image behind it had, so the bar is made here instead
    if (isTRUE(progress))
    {
        reporter <- newProgress(TRUE, length(index), list(perCall = 1, name = "voxels"))
        on.exit(reporter$close(), add = TRUE)
        progress <- function (done, total) reporter$report(done)
    }

    input <- maskedInput(x, selected, index, nLocations, dims, nSpatial)
    result <- imapply(input$data, input$margin, fun, ..., simplify = simplify,
                      threads = threads, progress = progress)

    scatterMasked(result, index, spatialDims, length(index), fill, x, nSpatial)
}

## A mask may be given as a logical array, a sparse image whose own mask is
## wanted, or anything numeric where non-zero means selected
asMaskVector <- function (mask, spatialDims)
{
    if (isSparseImage(mask))
        mask <- mask(mask)

    if (is.logical(mask))
        selected <- as.vector(mask)
    else if (is.numeric(mask))
        selected <- as.vector(mask) != 0
    else
        stop("Mask must be a logical array, a numeric array, or a sparse image")

    if (anyNA(selected))
        stop("Mask must not contain missing values")
    if (length(selected) != prod(spatialDims))
        stop("Mask must have one element per spatial location (", prod(spatialDims), ")")

    selected
}

## The selected values, and the margin to apply over. Both layouts hand the
## function the same vector; they differ only in which is cheaper to produce
maskedInput <- function (x, selected, index, nLocations, dims, nSpatial)
{
    elements <- if (nSpatial < length(dims)) prod(dims[-seq_len(nSpatial)]) else 1L

    ## The values a sparse image already holds, if they are the ones asked for
    if (isSparseImage(x) && identical(as.vector(mask(x)), selected))
        return(list(data = maskedMatrix(x), margin = 2L))

    values <- as.array(asDense(x))
    dim(values) <- c(nLocations, elements)
    list(data = values[index, , drop = FALSE], margin = 1L)
}

## Put the answers back where they came from, leaving fill everywhere else
scatterMasked <- function (result, index, spatialDims, nSelected, fill, x, nSpatial)
{
    if (!is.atomic(result) || length(result) %% nSelected != 0L)
        return(result)

    perLocation <- length(result) %/% nSelected
    full <- array(as.vector(fill, mode = typeof(result)), c(perLocation, prod(spatialDims)))
    full[, index] <- result

    if (perLocation == 1L)
    {
        dim(full) <- spatialDims
        if (isImage(x) && typeof(full) %in% c("logical", "integer", "double", "complex"))
            return(denseImage(full, template = x, spatial = nSpatial))
        return(full)
    }

    dim(full) <- c(perLocation, spatialDims)
    full
}
