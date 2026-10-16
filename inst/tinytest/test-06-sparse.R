## Sparse images.
##
## The contract is that a sparse image behaves exactly like the dense image it
## packs, so nearly everything here compares the two directly. What differs is
## how much is stored, and whether an operation can avoid unpacking.

set.seed(1)
dense <- array(0, c(8L, 8L, 4L))
dense[sample(256L, 60L)] <- rnorm(60L)
dense[3, 3, 1] <- NA                      # NA is not zero, so it must survive

other <- array(0, c(8L, 8L, 4L))
other[sample(256L, 40L)] <- rnorm(40L)

image <- denseImage(dense, voxelSize = c(2, 2, 2))
s <- asSparse(image)
t <- asSparse(denseImage(other))

## --- Construction and round trip -------------------------------------------

expect_true(isSparseImage(s))
expect_false(isDenseImage(s))
expect_equal(dim(s), c(8L, 8L, 4L))
expect_equal(length(s), 256L)
expect_equal(spatial(s), 3L)
expect_identical(geometry(s), geometry(image))

## Packing loses nothing at all, missing values included
expect_identical(as.array(s), dense)
expect_identical(as.array(asDense(s)), dense)
expect_true(is.na(as.array(s)[3, 3, 1]))

## Geometry travels with the image
expect_equal(voxelSize(s), c(2, 2, 2))
expect_equal(worldTransform(s), worldTransform(image))
expect_equal(voxelSize(asDense(s)), c(2, 2, 2))

## A location holding NA is stored, not folded away
expect_equal(sum(mask(s)), 61L)
expect_equal(sum(mask(s)), sum(apply(dense, 1:3, function (v) !identical(v, 0))))
expect_true(mask(s)[3, 3, 1])

## Sparseness is the proportion of locations holding nothing
expect_equal(sparseness(s), 1 - 61 / 256)

## asSparse on something already sparse is a no-op
expect_identical(asSparse(s), s)

## Every storage mode packs
for (mode in c("logical", "integer", "double", "complex"))
{
    a <- array(vector(mode, 64L), c(4L, 4L, 4L))
    a[c(2L, 7L, 30L)] <- as.vector(c(1, 1, 1), mode = mode)
    packed <- asSparse(denseImage(a))
    expect_identical(as.array(packed), a, info = mode)
    expect_identical(typeof(packed@values), mode, info = mode)
    expect_equal(sum(mask(packed)), 3L, info = mode)
}

## --- Values per location ---------------------------------------------------

## Sparsity is over locations: a location is present with all of its values, or
## absent with none of them
series <- array(0, c(4L, 4L, 4L, 5L))
series[2, 2, 2, ] <- rnorm(5)
series[3, 3, 3, 1] <- 1                    # only one value, but the location counts
packedSeries <- asSparse(denseImage(series))

expect_equal(sum(mask(packedSeries)), 2L)
expect_identical(as.array(packedSeries), series)
expect_equal(length(packedSeries@values), 2L * 5L)
expect_equal(dim(mask(packedSeries)), c(4L, 4L, 4L))

## --- Indexing --------------------------------------------------------------

## Lookup goes through the mask rather than materialising the image
expect_equal(s[1], dense[1])
expect_equal(s[c(1L, 50L, 200L)], dense[c(1L, 50L, 200L)])
expect_equal(s[], dense)

locations <- matrix(c(3L, 3L, 1L, 1L, 1L, 1L, 8L, 8L, 4L), ncol = 3L, byrow = TRUE)
expect_equal(s[locations], dense[locations])

## Full n-dimensional indexing agrees with the dense image
expect_equal(s[, , 1], dense[, , 1])
expect_equal(s[2:4, , 2], dense[2:4, , 2])

## --- Arithmetic ------------------------------------------------------------

## Each pair states the sparse expression and the dense one it must match, and
## whether the result should still be packed. An operation that leaves zero
## alone stays sparse; one that does not genuinely produces a dense image, and
## saying so beats re-packing at some arbitrary threshold
cases <- list(
    list(quote(s * 2),          quote(dense * 2),           TRUE),
    list(quote(s / 2),          quote(dense / 2),           TRUE),
    list(quote(2 * s),          quote(2 * dense),           TRUE),
    list(quote(0 - s),          quote(0 - dense),           TRUE),
    list(quote(s - s),          quote(dense - dense),       TRUE),
    list(quote(s * s),          quote(dense * dense),       TRUE),
    list(quote(s ^ 2),          quote(dense ^ 2),           TRUE),
    list(quote(s %% 2),         quote(dense %% 2),          TRUE),
    list(quote(s > 0),          quote(dense > 0),           TRUE),
    list(quote(s < 0),          quote(dense < 0),           TRUE),
    list(quote(abs(s)),         quote(abs(dense)),          TRUE),
    list(quote(sqrt(abs(s))),   quote(sqrt(abs(dense))),    TRUE),
    list(quote(sign(s)),        quote(sign(dense)),         TRUE),
    list(quote(round(s)),       quote(round(dense)),        TRUE),
    list(quote(s + 1),          quote(dense + 1),           FALSE),
    list(quote(1 + s),          quote(1 + dense),           FALSE),
    list(quote(s - 1),          quote(dense - 1),           FALSE),
    list(quote(exp(s)),         quote(exp(dense)),          FALSE),
    list(quote(cos(s)),         quote(cos(dense)),          FALSE),
    list(quote(s == 0),         quote(dense == 0),          FALSE),
    list(quote(s * t),          quote(dense * other),       TRUE),
    list(quote(s + t),          quote(dense + other),       TRUE),
    list(quote(s - t),          quote(dense - other),       TRUE),
    list(quote(s & t),          quote(dense & other),       TRUE),
    list(quote(s | t),          quote(dense | other),       TRUE),
    list(quote(s > t),          quote(dense > other),       TRUE),
    list(quote(s == t),         quote(dense == other),      FALSE))

for (case in cases)
{
    label <- deparse(case[[1L]])
    result <- eval(case[[1L]])
    expected <- eval(case[[2L]])

    expect_identical(as.array(result), expected, info = paste(label, "value"))
    expect_equal(isSparseImage(result), case[[3L]], info = paste(label, "stays packed"))
    expect_true(isSparseImage(result) || isDenseImage(result), info = paste(label, "is an image"))
}

## Geometry survives arithmetic
expect_equal(voxelSize(s * 2), c(2, 2, 2))
expect_equal(worldTransform(s * 2), worldTransform(image))
expect_equal(voxelSize(s + 1), c(2, 2, 2))

## Mixing a sparse and a dense image works in both directions
expect_identical(as.array(s * image), dense * dense)
expect_identical(as.array(image * s), dense * dense)

## Mismatched shapes are refused rather than silently recycled
expect_error(s * asSparse(denseImage(array(0, c(4L, 4L, 4L)))), "same dimensions")

## --- Tightening ------------------------------------------------------------

## An operation can turn stored values into zeros. The mask is re-derived
## afterwards, so sparseness() never reports a figure that has stopped being
## true. Only the NA location survives multiplication by zero
zeroed <- s * 0
expect_true(isSparseImage(zeroed))
expect_equal(sum(mask(zeroed)), 1L)
expect_true(is.na(as.array(zeroed)[3, 3, 1]))
expect_identical(as.array(zeroed), dense * 0)

## Multiplication annihilates where either operand is absent, so the stored
## locations collapse to the intersection without that being special-cased.
## NA is the exception, since NA * 0 is NA rather than zero, so the property
## is stated on images that hold none
cleanA <- asSparse(denseImage(replace(dense, is.na(dense), 0)))
cleanProduct <- cleanA * t
expect_true(all(mask(cleanProduct) <= (mask(cleanA) & mask(t))))
expect_true(sum(mask(cleanProduct)) < sum(mask(cleanA)))

## With an NA present, the location survives multiplication by an absent one,
## which is the arithmetic being right rather than the packing being wrong
product <- s * t
expect_true(mask(product)[3, 3, 1])
expect_false(mask(t)[3, 3, 1])
expect_true(all(mask(product) <= (mask(s) & mask(t)) | is.na(as.array(product))))

## --- Summaries -------------------------------------------------------------

## These are answered from the stored values plus one zero standing for every
## absent location, so they never unpack
expect_identical(sum(s, na.rm = TRUE), sum(dense, na.rm = TRUE))
expect_identical(max(s, na.rm = TRUE), max(dense, na.rm = TRUE))
expect_identical(min(s, na.rm = TRUE), min(dense, na.rm = TRUE))
expect_identical(range(s, na.rm = TRUE), range(dense, na.rm = TRUE))
expect_identical(prod(s, na.rm = TRUE), prod(dense, na.rm = TRUE))
expect_identical(any(s > 0, na.rm = TRUE), any(dense > 0, na.rm = TRUE))
expect_identical(all(s > 0, na.rm = TRUE), all(dense > 0, na.rm = TRUE))
expect_identical(mean(s, na.rm = TRUE), sum(dense, na.rm = TRUE) / length(dense))

## Missing values propagate as they would for the dense image
expect_identical(sum(s), sum(dense))
expect_true(is.na(sum(s)))

## A fully occupied image still summarises correctly, with no zero added
full <- asSparse(denseImage(array(rnorm(64), c(4L, 4L, 4L))))
expect_equal(sparseness(full), 0)
expect_identical(min(full), min(as.array(full)))
expect_identical(max(full), max(as.array(full)))

## --- Mask operations -------------------------------------------------------

expect_equal(imply:::maskCount(s@mask, 256), 61)
expect_identical(imply:::maskToLogical(imply:::maskFromLogical(c(TRUE, FALSE, TRUE)), 3), c(TRUE, FALSE, TRUE))
expect_error(imply:::maskFromLogical(c(TRUE, NA)), "must not contain missing")

union <- imply:::maskCombine(s@mask, t@mask, "union")
intersection <- imply:::maskCombine(s@mask, t@mask, "intersection")
expect_identical(imply:::maskToLogical(union, 256), as.vector(mask(s) | mask(t)))
expect_identical(imply:::maskToLogical(intersection, 256), as.vector(mask(s) & mask(t)))
expect_error(imply:::maskCombine(s@mask, t@mask, "nonsense"), "union")

## The bit count is right across word boundaries, which is where an off-by-one
## in the rank table would show up
for (n in c(1L, 63L, 64L, 65L, 127L, 128L, 129L, 1000L))
{
    present <- rep(c(TRUE, FALSE), length.out = n)
    packed <- imply:::maskFromLogical(present)
    expect_equal(imply:::maskCount(packed, n), sum(present), info = paste("count at", n))
    expect_identical(imply:::maskToLogical(packed, n), present, info = paste("round trip at", n))
}

## --- Degenerate cases ------------------------------------------------------

empty <- asSparse(denseImage(array(0, c(4L, 4L, 4L))))
expect_equal(sparseness(empty), 1)
expect_equal(length(empty@values), 0L)
expect_identical(as.array(empty), array(0, c(4L, 4L, 4L)))
expect_identical(sum(empty), 0)
expect_true(isSparseImage(empty * 2))

## A two-dimensional image
flat <- asSparse(denseImage(array(c(1, 0, 0, 2), c(2L, 2L)), spatial = 2L, voxelSize = c(1, 1)))
expect_equal(sum(mask(flat)), 2L)
expect_identical(as.array(flat), array(c(1, 0, 0, 2), c(2L, 2L)))

## --- Memory ----------------------------------------------------------------

## The point of packing: a mostly-empty image should cost markedly less than
## the dense one. A coordinate list would store several indices per value and
## could easily cost more
big <- array(0, c(64L, 64L, 64L))
big[sample(length(big), length(big) %/% 20L)] <- rnorm(length(big) %/% 20L)
packedBig <- asSparse(denseImage(big))

denseBytes <- as.numeric(object.size(big))
packedBytes <- as.numeric(object.size(packedBig@mask)) + as.numeric(object.size(packedBig@values))

expect_true(packedBytes < denseBytes / 3,
            info = sprintf("packed %.0f bytes against dense %.0f", packedBytes, denseBytes))
expect_identical(as.array(packedBig), big)
