#' Built-in reductions
#'
#' `imreduce()` computes a summary over the margins of an image without
#' calling back into R for each one. That matters for two reasons: there is no
#' per-call interpreter overhead, and (since nothing touches R) the loop can
#' run on worker threads. `imapply()` recognises the equivalent base functions
#' and routes calls here automatically where the answer has been shown to be
#' the same.
#'
#' Accumulation is done in double mode whatever the values are stored as, so
#' the answer does not depend on the storage type and error does not compound
#' with the number of values.
#'
#' Note that for the arithmetic reductions this may not agree with the base
#' equivalent to the last bit: R accumulates `sum()` and `mean()` in long
#' double, which on some platforms is wider than double. Differences are of
#' the order of `1e-16`. `imapply()` therefore routes automatically only where
#' the answers are provably identical (`min`, `max`, `range`, `which.min` and
#' `which.max`) and leaves the rest to be asked for deliberately.
#'
#' `"var"` is the variance of the values, as `var(as.vector(x))` would give.
#' It is not [stats::var()] applied to a sub-array, which for a matrix returns
#' a covariance matrix.
#'
#' @param x An image, or any array.
#' @param margin The dimensions to retain, as for [base::apply()].
#' @param what The reduction: one of `"sum"`, `"mean"`, `"min"`, `"max"`,
#'   `"range"`, `"prod"`, `"var"`, `"sd"`, `"which.min"`, `"which.max"`,
#'   `"any"`, `"all"` or `"countNA"`.
#' @param na.rm Whether to ignore missing values.
#' @param threads Number of threads, or `NULL` for the default. See
#'   [parallelism].
#' @return A vector or array, shaped over `margin` as [base::apply()] would
#'   shape it. `"range"` yields two values per call, so it gains a leading
#'   dimension.
#' @name imreduce
NULL

## Names as imreduce() knows them, against the functions imapply() will
## recognise as equivalent. Most are base, but var and sd come from stats
reductionSources <- c(sum = "base", mean = "base", min = "base", max = "base",
                      range = "base", prod = "base", var = "stats", sd = "stats",
                      which.min = "base", which.max = "base", any = "base", all = "base")

reductionNames <- names(reductionSources)

#' @rdname imreduce
#' @export
imreduce <- function (x, margin, what, na.rm = FALSE, threads = NULL)
{
    what <- match.arg(what, c(reductionNames, "countNA"))

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

    nThreads <- resolveThreads(threads)

    values <- if (isPackedImage(x))
        reduceOverMarginPacked(x@values, x@storageType, x@dims, margin, what,
                               x@slope, x@intercept, na.rm, nThreads, layoutArg(x))
    else if (isSparseImage(x))
        reduceOverMarginSparse(x@mask, x@values, x@dims, spatial(x), margin, what, na.rm, nThreads, layoutArg(x))
    else
        reduceOverMargin(unclassArray(x), margin, what, na.rm, nThreads)

    shapeReduction(values, dims[margin], margin, what, dimnames(x))
}

shapeReduction <- function (values, marginDims, margin, what, dimNames)
{
    marginNames <- if (!is.null(dimNames)) dimNames[margin] else NULL
    if (isAllNull(marginNames))
        marginNames <- NULL

    ## which.min and which.max report positions, which are whole numbers
    if (what %in% c("which.min", "which.max", "countNA"))
        values <- as.integer(values)
    else if (what %in% c("any", "all"))
        values <- as.logical(values)

    width <- if (what == "range") 2L else 1L

    if (width > 1L)
        return(array(values, c(width, marginDims),
                     if (is.null(marginNames)) NULL
                     else c(list(NULL), marginNames)))

    if (length(margin) == 1L)
    {
        if (!is.null(marginNames))
            names(values) <- marginNames[[1L]]
        return(values)
    }

    array(values, marginDims, marginNames)
}

## The base functions imapply() will route to imreduce(), keyed by the name
## imreduce() knows them under. Populated on load, since the functions have to
## be compared by identity
knownReductions <- new.env(parent = emptyenv())

registerReductions <- function ()
{
    for (name in reductionNames)
        assign(name, getExportedValue(reductionSources[[name]], name), envir = knownReductions)
    invisible(NULL)
}

## Only these are routed automatically, because only these are guaranteed to
## give bit-identical answers to base::apply(), which is a promise imapply()
## makes and should not quietly break. They involve comparison rather than
## arithmetic, so there is no accumulation order to differ over, and ties go
## to the first value exactly as which.min() and which.max() do.
##
## Deliberately absent:
##
##   sum, prod, mean, sd  arithmetic, and R accumulates in long double. On a
##                        platform where that is wider than double the answers
##                        differ in the last place; ours differ from mean()
##                        by around 1e-16 even here
##
##   var                  not the same function at all. var() of a matrix is a
##                        covariance matrix, whereas imreduce() gives the
##                        variance of the values, as var(as.vector(x)) would
##
##   any, all             base warns when coercing double data to logical, and
##                        routing would silently swallow the warning
##
## All of them remain available through imreduce() itself, which is documented
## as computing in double with its own accumulation order.
exactReductions <- c("min", "max", "range", "which.min", "which.max")

routableReduction <- function (fun, x, extras, simplify)
{
    if (!simplify)
        return(NULL)
    if (length(extras) > 0L && !identical(names(extras), "na.rm"))
        return(NULL)

    ## imreduce() always works in double, so routing integer or logical data
    ## would change the type of the answer: sum() of integers is an integer,
    ## and min() preserves its input type
    ## A packed image is always widened to double on the way out, whatever it
    ## is stored as, so it is routable regardless
    values <- if (isPackedImage(x)) 1.0 else if (isSparseImage(x)) x@values else x
    if (!is.double(values))
        return(NULL)

    for (name in exactReductions)
    {
        if (identical(fun, get(name, envir = knownReductions)))
            return(name)
    }

    NULL
}
