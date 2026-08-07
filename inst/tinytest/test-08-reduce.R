## Built-in reductions, and the routing that sends imapply() to them.

set.seed(1)
a <- array(rnorm(12L * 13L * 20L), c(12L, 13L, 20L))

## --- Agreement with the base equivalents ------------------------------------

pairs <- list(
    list("sum",       base::sum),
    list("mean",      base::mean),
    list("min",       base::min),
    list("max",       base::max),
    list("prod",      base::prod),
    list("sd",        stats::sd),
    list("which.min", base::which.min),
    list("which.max", base::which.max))

for (pair in pairs)
{
    for (margin in list(3L, c(1L, 2L)))
    {
        label <- paste(pair[[1L]], "margin", paste(margin, collapse = ","))
        expect_equal(as.vector(imreduce(a, margin, pair[[1L]])),
                     as.vector(apply(a, margin, pair[[2L]])), info = label)
    }
}

## range yields two values per call, so it gains a leading dimension
expect_equal(imreduce(a, 3, "range"), apply(a, 3, range))
expect_equal(dim(imreduce(a, c(1, 2), "range")), c(2L, 12L, 13L))
expect_equal(imreduce(a, c(1, 2), "range"), apply(a, c(1, 2), range))

## var is the variance of the values, which for a sub-array of more than one
## dimension is not what stats::var() computes
expect_equal(imreduce(a, 3, "var")[1L], var(as.vector(a[, , 1L])))
expect_equal(as.vector(imreduce(a, 3, "var")), as.vector(imreduce(a, 3, "sd"))^2)

## Logical reductions
logical3d <- array(rep(c(TRUE, FALSE), length.out = 240L), c(4L, 6L, 10L))
expect_equal(imreduce(logical3d, 3, "any"), apply(logical3d, 3, any))
expect_equal(imreduce(logical3d, 3, "all"), apply(logical3d, 3, all))
expect_true(is.logical(imreduce(logical3d, 3, "any")))
expect_true(is.integer(imreduce(a, 3, "which.max")))

## Positions are one-based, and ties go to the first value as R's do
ties <- array(c(1, 0, 0, 1), c(2L, 2L))
expect_equal(imreduce(ties, 2, "which.min"), apply(ties, 2, which.min))
expect_equal(imreduce(ties, 2, "which.max"), apply(ties, 2, which.max))

## --- Missing values --------------------------------------------------------

withNA <- a
withNA[1L, 1L, 1L] <- NA
withNA[2L, 2L, 5L] <- NA

for (what in c("sum", "mean", "min", "max", "prod", "sd"))
{
    fn <- if (what == "sd") stats::sd else get(what, baseenv())
    expect_equal(as.vector(imreduce(withNA, 3, what)),
                 as.vector(apply(withNA, 3, fn)), info = paste(what, "with NA"))
    expect_equal(as.vector(imreduce(withNA, 3, what, na.rm = TRUE)),
                 as.vector(apply(withNA, 3, function (v) fn(v, na.rm = TRUE))),
                 info = paste(what, "with NA, removed"))
}

expect_equal(imreduce(withNA, 3, "range"), apply(withNA, 3, range))
expect_equal(imreduce(withNA, 3, "range", na.rm = TRUE),
             apply(withNA, 3, range, na.rm = TRUE))

## countNA has no base equivalent, and counts what the others would discard
expect_equal(sum(imreduce(withNA, 3, "countNA")), 2L)
expect_equal(imreduce(withNA, 3, "countNA")[1L], 1L)

## --- Parallelism -----------------------------------------------------------

## Nothing here touches R, so this is the tier that can actually use threads.
## Dividing the work must not change the answer
for (what in c("sum", "mean", "min", "max", "sd", "which.max", "range"))
    expect_identical(imreduce(a, 3, what, threads = 2L), imreduce(a, 3, what, threads = 1L),
                     info = paste(what, "across threads"))

expect_identical(imreduce(a, c(1, 2), "sum", threads = 2L), imreduce(a, c(1, 2), "sum", threads = 1L))

## --- Other representations -------------------------------------------------

image <- denseImage(a)
packed <- asPacked(image, "float32")
sparse <- asSparse(denseImage(replace(a, abs(a) < 1, 0)))

expect_equal(imreduce(packed, 3, "sum"), imreduce(as.array(packed), 3, "sum"))
expect_equal(imreduce(packed, 3, "max"), imreduce(as.array(packed), 3, "max"))
expect_identical(imreduce(sparse, 3, "sum"), imreduce(as.array(sparse), 3, "sum"))
expect_identical(imreduce(sparse, 3, "which.max"), imreduce(as.array(sparse), 3, "which.max"))
expect_identical(imreduce(sparse, 3, "sum", threads = 2L), imreduce(sparse, 3, "sum", threads = 1L))

## --- Validation ------------------------------------------------------------

expect_error(imreduce(a, 3, "median"), "should be one of")
expect_error(imreduce(a, 4, "sum"), "out of range")
expect_error(imreduce(a, c(1, 1), "sum"), "repeated dimension")
expect_error(imreduce(a, integer(0), "sum"), "at least one dimension")

## --- Routing from imapply --------------------------------------------------

## Only reductions that are provably bit-identical to the base answer are
## routed. Comparison-based ones qualify; arithmetic ones do not, because R
## accumulates in long double and the last place can differ
for (what in c("min", "max", "range", "which.min", "which.max"))
{
    fn <- get(what, baseenv())
    expect_identical(imapply(a, 3, fn), apply(a, 3, fn), info = paste(what, "routed"))
    expect_identical(imapply(a, c(1, 2), fn), apply(a, c(1, 2), fn), info = paste(what, "routed, two margins"))
}

## Arithmetic reductions are left alone, so imapply stays exactly equal to
## base::apply rather than nearly equal
for (fn in list(base::sum, base::mean, base::prod, stats::sd))
    expect_identical(imapply(a, 3, fn), apply(a, 3, fn))

## The distinction is real: imreduce and base can differ in the last place
expect_equal(as.vector(imreduce(a, 3, "mean")), as.vector(apply(a, 3, mean)), tolerance = 1e-12)

## Routing respects na.rm, and declines anything else
expect_identical(imapply(withNA, 3, min, na.rm = TRUE), apply(withNA, 3, min, na.rm = TRUE))
expect_identical(imapply(withNA, 3, min), apply(withNA, 3, min))
expect_identical(imapply(a, 3, max, simplify = FALSE), apply(a, 3, max, simplify = FALSE))

## Integer and logical data are not routed, since imreduce works in double and
## min() preserves the type of its input
integers <- array(1:240, c(4L, 6L, 10L))
expect_identical(imapply(integers, 3, min), apply(integers, 3, min))
expect_true(is.integer(imapply(integers, 3, min)))
expect_identical(imapply(logical3d, 3, max), apply(logical3d, 3, max))

## A function that merely looks like a reduction is not routed
expect_identical(imapply(a, 3, function (v) min(v)), apply(a, 3, function (v) min(v)))
expect_identical(imapply(a, 3, function (v) min(v) + 1), apply(a, 3, function (v) min(v) + 1))

## Packed and sparse images route too, since both are read as double
expect_identical(imapply(packed, 3, max), imapply(as.array(packed), 3, max))
expect_identical(imapply(sparse, 3, max), imapply(as.array(sparse), 3, max))

## --- Speed -----------------------------------------------------------------

## The point of the compiled kernel: no interpreter call per sub-array
big <- array(rnorm(60L * 60L * 200L), c(60L, 60L, 200L))
compiled <- system.time(imreduce(big, 3, "sum"))[["elapsed"]]
interpreted <- system.time(apply(big, 3, sum))[["elapsed"]]
expect_true(compiled <= interpreted + 0.05,
            info = sprintf("imreduce %.3f s against apply %.3f s", compiled, interpreted))
