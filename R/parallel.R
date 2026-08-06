#' Parallelism
#'
#' Work is parallelised in one of two ways, according to what is being run.
#'
#' Compiled kernels run concurrently in process, using Grand Central Dispatch
#' where it is available and OpenMP otherwise. Which of these was compiled in
#' is reported by `parallelBackend()`.
#'
#' An R function cannot be called from a worker thread, because R is
#' single-threaded and its interpreter is not reentrant. Applying an R function
#' is therefore parallelised by forking the R session instead, dividing the
#' calls between workers. Forked workers share the image through copy-on-write,
#' so nothing large is duplicated. Forking is unavailable on Windows, where
#' such work runs serially; the alternative, a socket cluster, would copy the
#' whole image to every worker and so defeat the purpose.
#'
#' The number of threads may be given per call, or set globally with
#' `options(imply.threads = n)`. Leaving it unset lets compiled kernels use the
#' backend's own default while R functions run serially.
#'
#' @param threads A thread count, or `NULL` to consult `getOption("imply.threads")`.
#' @return `parallelBackend()` returns `"libdispatch"`, `"openmp"` or `"none"`.
#'   `canFork()` reports whether R functions can be run in parallel.
#' @name parallelism
NULL

#' @rdname parallelism
#' @export
parallelBackend <- function () parallelInfo()$backend

#' @rdname parallelism
#' @export
canFork <- function () .Platform$OS.type != "windows"

#' @rdname parallelism
#' @export
resolveThreads <- function (threads = NULL)
{
    if (is.null(threads))
        threads <- getOption("imply.threads")
    if (is.null(threads))
        return(0L)

    threads <- as.integer(threads)[1L]
    if (is.na(threads) || threads < 1L)
        stop("Number of threads must be a positive integer")
    threads
}

## Divide a number of calls into contiguous chunks, one per worker. The chunks
## are half-open and zero-based, matching what the compiled side expects
callChunks <- function (nCalls, nChunks)
{
    nChunks <- max(1L, min(as.integer(nChunks), nCalls))
    size <- ceiling(nCalls / nChunks)
    starts <- seq(0, nCalls - 1, by = size)
    lapply(starts, function (start) c(start, min(start + size, nCalls)))
}

## Reassemble what the workers produced. The pieces are in call order, so this
## is a concatenation; it only has to drop to a list if the workers disagreed
## about the shape of a result, which the serial path would also have done
combineParts <- function (parts)
{
    if (length(parts) == 1L)
        return(parts[[1L]])

    failed <- vapply(parts, inherits, NA, "try-error")
    if (any(failed))
        stop("Applying the function in parallel failed: ",
             conditionMessage(attr(parts[[which(failed)[1L]]], "condition")))

    flat <- !vapply(parts, function (p) p$isList, NA)
    lengths <- vapply(parts, function (p) p$elementLength, numeric(1))

    if (all(flat) && length(unique(lengths)) == 1L && !is.na(lengths[1L]))
        return(list(values = unlist(lapply(parts, `[[`, "values"), recursive = FALSE, use.names = FALSE),
                    elementLength = lengths[1L],
                    isList = FALSE))

    list(values = unlist(lapply(parts, partAsList), recursive = FALSE, use.names = FALSE),
         elementLength = NA_real_,
         isList = TRUE)
}

partAsList <- function (part)
{
    if (part$isList)
        return(part$values)

    n <- length(part$values) / part$elementLength
    if (n == 0)
        return(list())
    unname(split(part$values, rep(seq_len(n), each = part$elementLength)))
}

## Run the compiled loop over the whole call space, or over chunks of it in
## forked workers
runOverMargin <- function (x, margin, wrapped, callNames, simplify, nCalls, threads)
{
    if (threads <= 1L || nCalls <= 1L || !canFork())
        return(applyOverMargin(x, margin, wrapped, callNames, simplify))

    chunks <- callChunks(nCalls, threads)
    if (length(chunks) == 1L)
        return(applyOverMargin(x, margin, wrapped, callNames, simplify))

    parts <- parallel::mclapply(chunks, function (chunk)
        applyOverMargin(x, margin, wrapped, callNames, simplify, chunk[1L], chunk[2L]),
        mc.cores = length(chunks))

    combineParts(parts)
}
