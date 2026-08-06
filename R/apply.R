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
#' passes a slice's worth. This is forced for `voxelApply()`, since a location
#' on its own carries nothing to compute with, and the other two follow it so
#' that the three read as one progression.
#'
#' \preformatted{
#'   dim(x) = 4 x 5 x 6 x 10, spatial = 3
#'
#'   voxelApply(x, f)           f sees  10
#'   lineApply(x, f, axis = 1)  f sees  4 x 10
#'   sliceApply(x, f, axis = 3) f sees  4 x 5 x 10
#' }
#'
#' `axis` always names a spatial dimension, but means what the name of each
#' function implies: a line *runs along* its axis, whereas a slice is *cut
#' across* its axis. Use [imapply()] directly to iterate over a non-spatial
#' dimension, such as applying a function to each volume in a time series.
#'
#' @param x An image, or any array. No custom class is required.
#' @param margin The dimensions to retain, as for [base::apply()].
#' @param fun A function to apply.
#' @param ... Further arguments to `fun`.
#' @param simplify Whether to simplify the result to an array where possible.
#' @param threads Number of threads to use, or `NULL` to consult
#'   `getOption("imply.threads")`. See [parallelism].
#' @param axis For `lineApply()`, the axis lines run along; for `sliceApply()`,
#'   the axis slices are cut across.
#' @return For `imapply()`, as [base::apply()]. For `voxelApply()`, an image
#'   when the function returns a single value per location, otherwise an array.
#' @name imapply
NULL

#' @rdname imapply
#' @export
imapply <- function (x, margin, fun, ..., simplify = TRUE, threads = NULL)
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

    if (nCalls == 0L)
        return(applyToEmpty(x, wrapped, callDims, callNames, marginDims, marginNames))

    out <- runOverMargin(unclassArray(x), margin, wrapped, callNames, simplify,
                         nCalls, resolveThreads(threads))
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
## the empty result. This is what base::apply() does
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
voxelApply <- function (x, fun, ..., simplify = TRUE, threads = NULL)
{
    nSpatial <- spatial(x)
    dims <- dim(x)
    if (is.null(dims))
        dims <- length(x)

    if (nSpatial < 1L)
        stop("Image has no spatial dimensions to apply over")
    if (nSpatial == length(dims))
        stop("Image holds a single value at each location, so there is nothing to apply over")

    result <- imapply(x, seq_len(nSpatial), fun, ..., simplify = simplify, threads = threads)

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
lineApply <- function (x, fun, ..., axis = 1L, simplify = TRUE, threads = NULL)
{
    nSpatial <- spatial(x)
    axis <- checkAxis(axis, nSpatial)

    if (nSpatial < 2L)
        stop("An image with one spatial dimension is a single line, so apply the function to it directly")

    ## Every spatial dimension but the axis is retained, so fun sees one line
    ## running along the axis, together with the values at each of its
    ## locations. Lines never overlap, which is what makes this decomposition
    ## safe to parallelise
    imapply(x, seq_len(nSpatial)[-axis], fun, ..., simplify = simplify, threads = threads)
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
sliceApply <- function (x, fun, ..., axis = 3L, simplify = TRUE, threads = NULL)
{
    nSpatial <- spatial(x)

    if (nSpatial < 3L)
        stop("Slices need three spatial dimensions; use lineApply() for a two-dimensional image")
    axis <- checkAxis(axis, nSpatial)

    ## Only the axis is retained, so fun sees the plane cut across it, together
    ## with the values at each of its locations
    imapply(x, axis, fun, ..., simplify = simplify, threads = threads)
}
