#' Progress reporting
#'
#' The apply verbs are most useful when each call takes long enough to be
#' worth the machinery, which is exactly when it is worth knowing how far
#' through the operation is. Passing `progress = TRUE` draws a text bar;
#' passing a function of `(done, total)` reports however you like.
#'
#' Besides the percentage complete, the bar shows throughput in voxels per
#' second. Voxels are used rather than calls because they are comparable
#' across the verbs: a call is one location for [voxelApply()] but a whole
#' plane for [sliceApply()], whereas both sweep the same number of voxels. The
#' figure is a running average over the operation so far, rather than an
#' instantaneous rate that would flicker from one update to the next. A
#' function reporter is still given calls, since that is the unit its total is
#' in, and is free to convert.
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

## Returns NULL, or a pair of closures: one to report a count, one to tidy up.
## A function reporter is told about calls, since that is the unit it was
## given a total in; only the bar converts to voxels
newProgress <- function (progress, total, unit = NULL)
{
    if (is.null(progress) || isFALSE(progress))
        return(NULL)

    if (is.function(progress))
        return(list(report = function (done) progress(done, total),
                    close = function () invisible(NULL)))

    if (!isTRUE(progress))
        stop("progress must be TRUE, FALSE, or a function of (done, total)")

    if (is.null(unit))
        unit <- list(perCall = 1, name = "calls")
    started <- proc.time()[["elapsed"]]

    list(report = function (done) drawBar(done, total, unit, started),
         close = function () { drawBar(total, total, unit, started); cat("\n") })
}

## What a single call covers, for the rate. A call is handed every dimension
## outside the margin, so the locations it touches are the spatial ones among
## them: none for voxelApply(), a line's worth for lineApply(), and so on
progressUnit <- function (x, dims, margin)
{
    nSpatial <- spatial(x)
    if (nSpatial < 1L)
        return(list(perCall = 1, name = "calls"))
    list(perCall = prod(dims[setdiff(seq_len(nSpatial), margin)]), name = "voxels")
}

## Drawn by hand rather than with utils::txtProgressBar(), whose label argument
## is ignored, so there would be nowhere to put the rate. Every line is the
## same width, so a carriage return is enough to overwrite the last one
drawBar <- function (done, total, unit, started)
{
    fraction <- if (total > 0) min(1, done / total) else 1
    elapsed <- proc.time()[["elapsed"]] - started
    rate <- if (elapsed > 0) done * unit$perCall / elapsed else NA_real_

    suffix <- sprintf("%3.0f%%  %s %s/s", 100 * fraction, formatRate(rate), unit$name)
    width <- max(20L, min(as.integer(getOption("width", 80L)), 100L))
    barWidth <- max(10L, width - nchar(suffix) - 6L)
    filled <- as.integer(round(fraction * barWidth))

    cat("\r  |", strrep("=", filled), strrep(" ", barWidth - filled), "| ", suffix, sep = "")
    utils::flush.console()
    invisible(NULL)
}

## Three significant figures at most, so the width of the line never changes
## enough to leave debris behind
formatRate <- function (rate)
{
    if (!is.finite(rate))
        return("    --")
    scale <- min(4L, sum(rate >= c(1e3, 1e6, 1e9, 1e12)))
    scaled <- rate / 1000^scale
    sprintf("%6s", paste0(if (scaled >= 100) sprintf("%.0f", scaled) else sprintf("%.1f", scaled),
                          c("", "k", "M", "G", "T")[scale + 1L]))
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
