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

## --- Interruption ----------------------------------------------------------

## Rcpp::checkUserInterrupt() runs R_CheckUserInterrupt() inside
## R_ToplevelExec, so R's longjmp is contained and an ordinary C++ exception
## is thrown instead, letting the sinks and walkers be destroyed on the way
## out. Raising a real SIGINT would kill the test process, so what is checked
## here is that adding the checks changed no answer, and that the slabbed
## reduction path agrees with the unslabbed one.

## The apply loop looks for an interrupt every hundred calls, which has to
## divide the work without disturbing it
manyCalls <- array(rnorm(4L * 4L * 250L), c(4L, 4L, 250L))
expect_identical(imapply(manyCalls, 3, function (v) sum(v)), apply(manyCalls, 3, sum))
expect_identical(imapply(manyCalls, c(1, 2), function (v) mean(v)), apply(manyCalls, c(1, 2), mean))

## A reduction large enough to be split into slabs, so that it can be
## interrupted between them, must give the same answer as a small one that is
## not split. The threshold is on total elements, so this crosses it
large <- array(rnorm(120L * 120L * 800L), c(120L, 120L, 800L))
expect_true(length(large) > 1e7)
expect_equal(as.vector(imreduce(large, 3, "sum")), as.vector(apply(large, 3, sum)))
expect_equal(as.vector(imreduce(large, 3, "max")), as.vector(apply(large, 3, max)))
expect_identical(imreduce(large, 3, "sum", threads = 2L), imreduce(large, 3, "sum", threads = 1L))

## ...and slabbing must not disturb a reduction whose result is wider than one
## value per call, where the output positions matter
expect_equal(imreduce(large, 3, "range"), apply(large, 3, range))
expect_equal(imreduce(large, 3, "which.max"), apply(large, 3, which.max))

## Below the threshold the work is done in one go, and still agrees
small <- array(rnorm(20L * 20L * 50L), c(20L, 20L, 50L))
expect_true(length(small) < 1e7)
expect_equal(as.vector(imreduce(small, 3, "sum")), as.vector(apply(small, 3, sum)))
