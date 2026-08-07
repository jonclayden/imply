#' Progress reporting
#'
#' The apply verbs are most useful when each call takes long enough to be
#' worth the machinery, which is exactly when it is worth knowing how far
#' through the operation is. Passing `progress = TRUE` draws a text bar;
#' passing a function of `(done, total)` reports however you like.
#'
#' Under forked parallelism the work is divided into batches, and the bar
#' advances as each batch completes rather than as each call does. A worker
#' cannot report on its own behalf: it is a separate process, and several of
#' them writing to one console would interleave. The batching costs one extra
#' fork per batch, which is negligible against work slow enough to want a
#' progress bar in the first place.
#'
#' @param progress `FALSE` for none, `TRUE` for a text bar, or a function
#'   called with the number of calls completed and the total.
#' @param total The number of calls that will be made.
#' @name progress
NULL

## Returns NULL, or a pair of closures: one to report a count, one to tidy up
newProgress <- function (progress, total)
{
    if (is.null(progress) || isFALSE(progress))
        return(NULL)

    if (is.function(progress))
        return(list(report = function (done) progress(done, total),
                    close = function () invisible(NULL)))

    if (!isTRUE(progress))
        stop("progress must be TRUE, FALSE, or a function of (done, total)")

    bar <- utils::txtProgressBar(min = 0, max = max(total, 1), style = 3)
    list(report = function (done) utils::setTxtProgressBar(bar, done),
         close = function () { utils::setTxtProgressBar(bar, total); close(bar); cat("\n") })
}

## How often the compiled loop should call back. Around a hundred updates is
## smooth without the callback itself becoming part of the cost
reportInterval <- function (total) max(1L, as.integer(total %/% 100L))

## Contiguous zero-based half-open ranges covering [from, to)
rangeChunks <- function (from, to, count)
{
    n <- to - from
    count <- max(1L, min(as.integer(count), n))
    size <- ceiling(n / count)
    starts <- seq(from, to - 1, by = size)
    lapply(starts, function (start) c(start, min(start + size, to)))
}
