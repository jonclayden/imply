## Progress reporting.
##
## The substantive requirement is that reporting never changes the answer, and
## that the counts it reports are monotonic, absolute, and end at the total.

set.seed(1)
x <- array(rnorm(10L * 10L * 40L), c(10L, 10L, 40L))

## A reporter that simply records what it was told
recorder <- function ()
{
    seen <- numeric(0)
    totals <- numeric(0)
    list(fn = function (done, total) {
             seen <<- c(seen, done)
             totals <<- c(totals, total)
         },
         seen = function () seen,
         totals = function () totals)
}

## --- Off by default --------------------------------------------------------

expect_silent(imapply(x, 3, sum))
expect_identical(imapply(x, 3, function (v) mean(v)),
                 imapply(x, 3, function (v) mean(v), progress = FALSE))

## --- A function reporter ---------------------------------------------------

r <- recorder()
plain <- imapply(x, 3, function (v) mean(v))
reported <- imapply(x, 3, function (v) mean(v), progress = r$fn)

## Reporting must not perturb the result
expect_identical(reported, plain)

## Counts rise, never repeat backwards, and finish at the number of calls
expect_true(length(r$seen()) > 1L)
expect_false(is.unsorted(r$seen()))
expect_equal(tail(r$seen(), 1L), 40)
expect_true(all(r$seen() >= 1 & r$seen() <= 40))

## The total is passed through unchanged every time
expect_equal(unique(r$totals()), 40)

## --- The text bar ----------------------------------------------------------

## Drawn on the console, so it is the side effect rather than the value that
## is checked
output <- capture.output(result <- imapply(x, 3, function (v) mean(v), progress = TRUE),
                         type = "output")
expect_identical(result, plain)
expect_true(any(nzchar(output)))
expect_true(any(grepl("100%", output, fixed = TRUE)))

expect_error(imapply(x, 3, sum, progress = "yes"), "TRUE, FALSE, or a function")
expect_error(imapply(x, 3, sum, progress = 1:3), "TRUE, FALSE, or a function")

## --- Under forked parallelism ----------------------------------------------

## A worker cannot report on its own behalf, so the work is batched and the
## parent reports between batches. The counts must still be absolute and
## monotonic, and the answer must be the serial one
if (canFork())
{
    p <- recorder()
    inParallel <- imapply(x, 3, function (v) mean(v), threads = 2L, progress = p$fn)

    expect_identical(inParallel, plain)
    expect_true(length(p$seen()) > 1L)
    expect_false(is.unsorted(p$seen()))
    expect_equal(tail(p$seen(), 1L), 40)
    expect_equal(unique(p$totals()), 40)

    ## Batching must not change how the results are assembled, including for
    ## the awkward shapes the workers might disagree about
    for (fn in list(function (v) range(v),
                    function (v) seq_len(1L + (round(abs(sum(v))) %% 3L)),
                    function (v) NULL))
    {
        expect_identical(imapply(x, 3, fn, threads = 2L, progress = function (d, t) NULL),
                         imapply(x, 3, fn, threads = 1L),
                         info = "batched parallel with progress matches serial")
    }
}

## --- The wrappers ----------------------------------------------------------

image <- denseImage(array(rnorm(6L * 6L * 6L * 10L), c(6L, 6L, 6L, 10L)))

for (call in list(quote(voxelApply(image, function (v) mean(v), progress = counter)),
                  quote(lineApply(image, function (v) sum(v), axis = 1, progress = counter)),
                  quote(sliceApply(image, function (v) sum(v), axis = 3, progress = counter))))
{
    calls <- 0L
    counter <- function (done, total) calls <<- calls + 1L
    invisible(eval(call))
    expect_true(calls > 0L, info = deparse(call))
}

## The result is the same with and without
expect_identical(voxelApply(image, function (v) mean(v), progress = function (d, t) NULL),
                 voxelApply(image, function (v) mean(v)))

## --- Other representations -------------------------------------------------

packed <- asPacked(image, "float32")
sparse <- asSparse(denseImage(replace(as.array(image), abs(as.array(image)) < 1, 0)))

for (compact in list(packed, sparse))
{
    calls <- 0L
    quiet <- imapply(compact, 4, function (v) mean(v))
    noisy <- imapply(compact, 4, function (v) mean(v),
                     progress = function (done, total) calls <<- calls + 1L)
    expect_identical(noisy, quiet)
    expect_true(calls > 0L)
}

## --- Interval --------------------------------------------------------------

## Roughly a hundred updates, so the callback is never a meaningful part of
## the cost, and never fewer than one call apart
expect_equal(imply:::reportInterval(10000), 100L)
expect_equal(imply:::reportInterval(50), 1L)
expect_equal(imply:::reportInterval(1), 1L)

## Ranges tile the call space exactly once, in order
for (n in c(1L, 7L, 100L))
{
    for (k in c(1L, 3L, 50L))
    {
        pieces <- imply:::rangeChunks(0, n, k)
        starts <- vapply(pieces, `[`, numeric(1), 1L)
        ends <- vapply(pieces, `[`, numeric(1), 2L)
        expect_equal(starts[1L], 0, info = paste(n, k))
        expect_equal(ends[length(ends)], n, info = paste(n, k))
        expect_equal(starts[-1L], head(ends, -1L), info = paste(n, k))
    }
}

## An offset range works the same way, which is what the batched parallel path
## relies on
pieces <- imply:::rangeChunks(10, 25, 4)
expect_equal(vapply(pieces, `[`, numeric(1), 1L)[1L], 10)
expect_equal(tail(vapply(pieces, `[`, numeric(1), 2L), 1L), 25)
