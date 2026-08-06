## imapply() and its wrappers.
##
## The contract is that imapply() gives exactly what base::apply() gives, so
## most of these tests are equivalence checks rather than fixed expectations.
## What differs is how much memory is used getting there.

named <- function (a) {
    dimnames(a) <- lapply(dim(a), function (k) paste0("d", seq_len(k)))
    a
}

set.seed(1)
arrays <- list(
    real       = array(rnorm(120), c(4L, 5L, 6L)),
    integer    = array(1:120, c(4L, 5L, 6L)),
    logical    = array(rep(c(TRUE, FALSE), 60), c(4L, 5L, 6L)),
    complex    = array(complex(real = 1:24, imaginary = 24:1), c(2L, 3L, 4L)),
    fourD      = array(rnorm(240), c(4L, 5L, 6L, 2L)),
    vector     = array(rnorm(60), 60L),
    matrix     = array(rnorm(120), c(4L, 30L)),
    withNames  = named(array(rnorm(120), c(4L, 5L, 6L))),
    zeroExtent = array(numeric(0), c(4L, 0L, 6L)))

## Functions chosen to exercise each branch of the result shaping: scalar,
## fixed-length vector, ragged, NULL, non-atomic, named, and attribute-carrying
functions <- list(
    sum       = sum,
    mean      = mean,
    range     = range,
    identity  = function (v) v,
    length    = length,
    firstTwo  = function (v) v[1:2],
    null      = function (v) NULL,
    character = function (v) "a",
    ragged    = function (v) seq_len(1L + (length(v) %% 3L)),
    named     = function (v) c(a = 1, b = 2),
    namedVary = function (v) setNames(1:2, c(letters[1L + (length(v) %% 5L)], "b")),
    list      = function (v) list(v),
    matrix    = function (v) matrix(1:4, 2L),
    quantile  = function (v) quantile(as.numeric(v)),
    usesDim   = function (v) dim(v),
    usesNames = function (v) names(v)[1L])

comparisons <- 0L
for (a in arrays)
{
    nDims <- length(dim(a))
    margins <- c(as.list(seq_len(nDims)),
                 if (nDims > 1L) list(c(1L, 2L)),
                 if (nDims > 2L) list(c(1L, 3L), c(3L, 1L), c(2L, 3L)))

    for (name in names(functions))
    {
        for (margin in margins)
        {
            for (simplify in c(TRUE, FALSE))
            {
                label <- paste0(name, " margin=", paste(margin, collapse = ","),
                                " simplify=", simplify, " dim=", paste(dim(a), collapse = "x"))
                reference <- suppressWarnings(tryCatch(
                    apply(a, margin, functions[[name]], simplify = simplify),
                    error = function (e) paste("error:", conditionMessage(e))))
                result <- suppressWarnings(tryCatch(
                    imapply(a, margin, functions[[name]], simplify = simplify),
                    error = function (e) paste("error:", conditionMessage(e))))

                expect_identical(result, reference, info = label)
                comparisons <- comparisons + 1L
            }
        }
    }
}

## Guard against the loop silently doing nothing
expect_true(comparisons > 1000L)

## --- Extra arguments -------------------------------------------------------

withNA <- array(c(NA, rnorm(119)), c(4L, 5L, 6L))
expect_identical(imapply(withNA, 3, mean, na.rm = TRUE), apply(withNA, 3, mean, na.rm = TRUE))
expect_identical(imapply(withNA, 3, mean), apply(withNA, 3, mean))
expect_identical(imapply(withNA, c(1, 2), quantile, probs = 0.5, na.rm = TRUE),
                 apply(withNA, c(1, 2), quantile, probs = 0.5, na.rm = TRUE))

## A function given by name, as base::apply accepts
expect_identical(imapply(arrays$real, 3, "sum"), apply(arrays$real, 3, "sum"))

## --- Validation ------------------------------------------------------------

x <- arrays$real
expect_error(imapply(x, 4, sum), "out of range")
expect_error(imapply(x, 0, sum), "out of range")
expect_error(imapply(x, c(1, 1), sum), "repeated dimension")
expect_error(imapply(x, integer(0), sum), "at least one dimension")
expect_error(imapply(x, NA_integer_, sum), "out of range")

## --- voxelApply ------------------------------------------------------------

image <- denseImage(array(rnorm(4 * 5 * 6 * 10), c(4L, 5L, 6L, 10L)), pixdim = c(2, 2, 3))

## Applying over the values at each location is applying over the spatial
## margins, and a single value per location is itself an image
means <- voxelApply(image, mean)
expect_true(isDenseImage(means))
expect_equal(dim(means), c(4L, 5L, 6L))
expect_equal(pixdim(means), c(2, 2, 3))
expect_equal(xform(means), xform(image))
expect_equal(as.array(means), apply(as.array(image), 1:3, mean))

## More than one value per location is an array, not an image, because the
## result no longer has the shape of the space
ranges <- voxelApply(image, range)
expect_false(isDenseImage(ranges))
expect_equal(dim(ranges), c(2L, 4L, 5L, 6L))
expect_equal(ranges, apply(as.array(image), 1:3, range))

## The function really does see the whole series at each location
expect_equal(unique(as.vector(voxelApply(image, length))), 10L)

## A two-dimensional image with a series at each location
slicesOverTime <- denseImage(array(rnorm(4 * 5 * 7), c(4L, 5L, 7L)), spatial = 2L, pixdim = c(1, 1))
expect_equal(dim(voxelApply(slicesOverTime, mean)), c(4L, 5L))
expect_equal(unique(as.vector(voxelApply(slicesOverTime, length))), 7L)

expect_error(voxelApply(denseImage(array(0, c(4L, 5L, 6L))), mean), "nothing to apply over")

## --- lineApply -------------------------------------------------------------

## The spatial unit grows from a location to a line to a slice, and the values
## held at each location travel with it. So on this image, which has ten
## values per location, a line along axis 1 is handed over as 4 x 10
plain <- as.array(image)
for (axis in 1:3)
{
    expect_identical(lineApply(image, sum, axis = axis),
                     apply(plain, seq_len(3)[-axis], sum),
                     info = paste("lineApply along axis", axis))
    expect_equal(dim(lineApply(image, sum, axis = axis)), dim(image)[seq_len(3)[-axis]])
    expect_equal(lineApply(image, dim, axis = axis)[, 1, 1],
                 c(dim(image)[axis], 10L),
                 info = paste("unit shape along axis", axis))
}

## Consistency of the progression: each verb differs from the next only in how
## much of the space it hands over
expect_equal(dim(voxelApply(image, identity)), c(10L, 4L, 5L, 6L))
expect_equal(lineApply(image, dim, axis = 1)[, 1, 1], c(4L, 10L))
expect_equal(sliceApply(image, dim, axis = 3)[, 1], c(4L, 5L, 10L))

## The two-dimensional case, where a slice is degenerate but a line is exactly
## what is wanted. With no values per location, a line is a bare vector
flat <- denseImage(array(rnorm(4 * 5), c(4L, 5L)), spatial = 2L, pixdim = c(1, 1))
expect_equal(length(lineApply(flat, sum, axis = 1)), 5L)
expect_equal(length(lineApply(flat, sum, axis = 2)), 4L)
expect_equal(unique(as.vector(lineApply(flat, length, axis = 1))), 4L)
expect_equal(unique(as.vector(lineApply(flat, length, axis = 2))), 5L)
expect_identical(lineApply(flat, sum, axis = 1), apply(as.array(flat), 2, sum))
expect_null(lineApply(flat, dim, axis = 1))

## A two-dimensional image with a series at each location: the line brings its
## series along
flatSeries <- denseImage(array(rnorm(4 * 5 * 7), c(4L, 5L, 7L)), spatial = 2L, pixdim = c(1, 1))
expect_equal(lineApply(flatSeries, dim, axis = 1)[, 1], c(4L, 7L))
expect_identical(lineApply(flatSeries, sum, axis = 1), apply(as.array(flatSeries), 2, sum))

## Results longer than one value stack on the leading dimension
expect_equal(dim(lineApply(image, range, axis = 1)), c(2L, 5L, 6L))
expect_identical(lineApply(image, range, axis = 1), apply(plain, c(2, 3), range))

## Extra arguments and simplify are forwarded
withNA <- plain
withNA[1, 1, 1, 1] <- NA
expect_identical(lineApply(withNA, mean, axis = 1, na.rm = TRUE),
                 apply(withNA, c(2, 3), mean, na.rm = TRUE))
expect_identical(lineApply(withNA, sum, axis = 1, simplify = FALSE),
                 apply(withNA, c(2, 3), sum, simplify = FALSE))

## Plain arrays behave identically
expect_identical(lineApply(plain, sum, axis = 2), lineApply(image, sum, axis = 2))

## The axis must be a spatial one; iterating over time is imapply()'s job
expect_error(lineApply(image, sum, axis = 4), "spatial dimension")
expect_error(lineApply(image, sum, axis = 9), "spatial dimension")
expect_error(lineApply(image, sum, axis = c(1, 2)), "spatial dimension")
expect_error(lineApply(image, sum, axis = NA), "spatial dimension")
expect_error(lineApply(denseImage(rnorm(10)), sum, axis = 1), "single line")

## --- sliceApply ------------------------------------------------------------

expect_equal(sliceApply(image, sum, axis = 3), apply(as.array(image), 3, sum))
expect_equal(sliceApply(image, sum, axis = 1), apply(as.array(image), 1, sum))
expect_equal(sliceApply(image, sum), apply(as.array(image), 3, sum))

## Each slice arrives with the remaining dimensions intact
expect_equal(sliceApply(image, dim, axis = 1)[, 1], c(5L, 6L, 10L))

expect_error(sliceApply(image, sum, axis = 4), "spatial dimension")
expect_error(sliceApply(image, sum, axis = 9), "spatial dimension")
expect_error(sliceApply(image, sum, axis = c(1, 2)), "spatial dimension")

## Slices need three spatial dimensions, which is what lineApply() is for
expect_error(sliceApply(flat, sum), "lineApply")

## --- Plain arrays and images agree -----------------------------------------

## No custom class is ever required, and using one changes nothing
expect_identical(voxelApply(as.array(image), mean), apply(as.array(image), 1:3, mean))
expect_identical(imapply(image, 4, sum), imapply(as.array(image), 4, sum))
expect_identical(sliceApply(as.array(image), sum, axis = 3), sliceApply(image, sum, axis = 3))

## --- Memory ----------------------------------------------------------------

## The point of the exercise: base::apply() permutes the whole array into a
## fresh copy before looping, so it allocates about twice the input. imapply()
## gathers through the stride vector instead and allocates about once
big <- array(rnorm(2e6), c(100L, 100L, 200L))
inputSize <- as.numeric(object.size(big))

measure <- function (expr) {
    invisible(gc(full = TRUE))
    before <- gc(reset = TRUE, full = TRUE)
    force(expr)
    after <- gc(full = TRUE)
    (after["Vcells", "max used"] - before["Vcells", "used"]) * 8
}

baseCost <- measure(apply(big, 3, mean))
implyCost <- measure(imapply(big, 3, mean))

expect_true(baseCost > 1.5 * inputSize,
            info = sprintf("base::apply allocated %.2f x input", baseCost / inputSize))
expect_true(implyCost < 1.25 * inputSize,
            info = sprintf("imapply allocated %.2f x input", implyCost / inputSize))
expect_true(implyCost < 0.75 * baseCost)
