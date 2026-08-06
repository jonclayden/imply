## Parallelism, both tiers.
##
## Tests stay at two threads, which is what CRAN permits during checks. The
## substantive requirement is that dividing the work never changes the answer.

## --- Backend reporting -----------------------------------------------------

info <- imply:::parallelInfo()
expect_true(info$backend %in% c("libdispatch", "openmp", "none"))
expect_identical(parallelBackend(), info$backend)
expect_identical(info$available, info$backend != "none")
expect_true(is.logical(canFork()) && length(canFork()) == 1L)

## --- Thread resolution -----------------------------------------------------

expect_equal(resolveThreads(NULL), 0L)
expect_equal(resolveThreads(4), 4L)
expect_equal(resolveThreads(1), 1L)
expect_error(resolveThreads(0), "positive integer")
expect_error(resolveThreads(-2), "positive integer")
expect_error(resolveThreads(NA), "positive integer")

## The global option is consulted when nothing is passed
old <- getOption("imply.threads")
options(imply.threads = 3L)
expect_equal(resolveThreads(NULL), 3L)
expect_equal(resolveThreads(2), 2L)          # an explicit value still wins
options(imply.threads = old)

## --- Work division ---------------------------------------------------------

## Chunks are what the backend iterates over, which is how a requested thread
## count comes to mean something under libdispatch. It offers no width control
## of its own, so the bound has to come from there being only so many chunks
if (info$available)
{
    expect_equal(imply:::chunkPartition(1000, 1), 1)
    expect_equal(imply:::chunkPartition(1000, 2), 2)
    expect_equal(imply:::chunkPartition(1000, 8), 8)

    ## Never more chunks than there are items to divide
    expect_equal(imply:::chunkPartition(3, 8), 3)
    expect_equal(imply:::chunkPartition(1, 8), 1)
}
expect_equal(imply:::chunkPartition(0, 4), 0)

## The R-side division covers the call space exactly once, in order
for (nCalls in c(1L, 2L, 7L, 100L))
{
    for (nChunks in c(1L, 2L, 3L, 8L, 50L))
    {
        chunks <- imply:::callChunks(nCalls, nChunks)
        starts <- vapply(chunks, `[`, numeric(1), 1L)
        ends <- vapply(chunks, `[`, numeric(1), 2L)

        expect_equal(starts[1L], 0, info = paste(nCalls, nChunks))
        expect_equal(ends[length(ends)], nCalls, info = paste(nCalls, nChunks))
        expect_equal(starts[-1L], head(ends, -1L), info = paste(nCalls, nChunks))
        expect_true(length(chunks) <= max(1L, min(nChunks, nCalls)), info = paste(nCalls, nChunks))
    }
}

## --- Compiled kernels ------------------------------------------------------

## Each worker owns a disjoint range of the output, so the result must not
## depend on how many workers there were
set.seed(1)
x <- array(rnorm(8 * 9 * 10 * 3), c(8L, 9L, 10L, 3L))
reference <- aperm(x, c(3, 1, 4, 2))

for (threads in c(1L, 2L))
{
    expect_identical(imply:::permuteView(x, c(3, 1, 4, 2), threads = threads), reference,
                     info = paste("permuteView with", threads, "threads"))
}

## Including on the runtime-dimensionality path
expect_identical(imply:::permuteView(x, c(3, 1, 4, 2), forceDynamic = TRUE, threads = 2L), reference)

## --- R callbacks -----------------------------------------------------------

## An R function cannot run on a worker thread, so this tier forks instead.
## Results must still be identical to the serial ones, whatever the shape
set.seed(2)
y <- array(rnorm(12 * 13 * 20), c(12L, 13L, 20L))

functions <- list(
    scalar    = mean,
    vector    = range,
    ragged    = function (v) seq_len(1L + (round(abs(sum(v))) %% 3L)),
    mixedType = function (v) if (mean(v) > 0) 1L else 2.5,
    null      = function (v) NULL,
    named     = function (v) c(low = min(v), high = max(v)),
    character = function (v) as.character(round(mean(v), 3)))

for (name in names(functions))
{
    for (margin in list(3L, c(1L, 2L)))
    {
        label <- paste(name, "margin", paste(margin, collapse = ","))
        serial <- imapply(y, margin, functions[[name]], threads = 1L)

        expect_identical(imapply(y, margin, functions[[name]], threads = 2L), serial,
                         info = paste(label, "parallel matches serial"))
        expect_identical(serial, apply(y, margin, functions[[name]]),
                         info = paste(label, "matches base::apply"))
    }
}

## simplify = FALSE, which forces the list path in every worker
expect_identical(imapply(y, 3, mean, threads = 2L, simplify = FALSE),
                 apply(y, 3, mean, simplify = FALSE))

## Extra arguments survive the trip to a worker
withNA <- y
withNA[1, 1, 1] <- NA
expect_identical(imapply(withNA, 3, mean, na.rm = TRUE, threads = 2L),
                 apply(withNA, 3, mean, na.rm = TRUE))

## The wrappers pass the thread count through
image <- denseImage(array(rnorm(6 * 7 * 8 * 5), c(6L, 7L, 8L, 5L)))
expect_identical(voxelApply(image, mean, threads = 2L), voxelApply(image, mean, threads = 1L))
expect_identical(lineApply(image, sum, axis = 1, threads = 2L), lineApply(image, sum, axis = 1, threads = 1L))
expect_identical(sliceApply(image, sum, axis = 3, threads = 2L), sliceApply(image, sum, axis = 3, threads = 1L))
expect_true(isDenseImage(voxelApply(image, mean, threads = 2L)))

## The option applies to the apply verbs too
options(imply.threads = 2L)
expect_identical(imapply(y, 3, mean), apply(y, 3, mean))
options(imply.threads = old)

## --- Partial ranges --------------------------------------------------------

## The compiled loop can be asked for part of the call space, which is what
## the workers are given. Concatenating the parts must reproduce the whole
whole <- imply:::applyOverMargin(y, 3L, function (v) mean(v))
first <- imply:::applyOverMargin(y, 3L, function (v) mean(v), from = 0, to = 8)
second <- imply:::applyOverMargin(y, 3L, function (v) mean(v), from = 8, to = 20)

expect_equal(length(first$values), 8L)
expect_equal(length(second$values), 12L)
expect_equal(c(first$values, second$values), whole$values)

## An empty range yields nothing rather than failing
empty <- imply:::applyOverMargin(y, 3L, function (v) mean(v), from = 5, to = 5)
expect_equal(length(empty$values), 0L)
