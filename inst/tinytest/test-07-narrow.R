## Narrow storage types.
##
## The contract is that a packed image behaves like the double image it
## approximates, to within the precision of the type it is stored in, while
## occupying a known fraction of the memory.

set.seed(1)
values <- array(rnorm(8L * 8L * 4L), c(8L, 8L, 4L))
image <- denseImage(values, voxelSize = c(2, 2, 2))

## --- float32 ---------------------------------------------------------------

f <- asPacked(image, "float32")

expect_true(isPackedImage(f))
expect_false(isDenseImage(f))
expect_equal(dim(f), c(8L, 8L, 4L))
expect_equal(length(f), 256L)
expect_equal(storageType(f), "float32")
expect_equal(f@slope, 1)
expect_equal(f@intercept, 0)

## Exactly half the memory of a double image, which is the point
expect_equal(length(f@values), 256L * 4L)
expect_equal(length(f@values), length(values) * 8L / 2L)

## Values agree to single precision
expect_true(max(abs(as.array(f) - values)) < 1e-6)
expect_equal(as.array(f), values, tolerance = 1e-6)

## Geometry travels with the image
expect_equal(voxelSize(f), c(2, 2, 2))
expect_equal(worldTransform(f), worldTransform(image))
expect_equal(spatial(f), 3L)

## --- Integer types and calibration -----------------------------------------

## Whole numbers already inside the range are stored as they are, so the
## packing is exact and the stored values stay readable
counts <- denseImage(array(sample(0:1000, 256L, TRUE), c(8L, 8L, 4L)))
i16 <- asPacked(counts, "int16")

expect_equal(i16@slope, 1)
expect_equal(i16@intercept, 0)
expect_equal(as.array(i16), as.array(counts) * 1.0)
expect_equal(length(i16@values), 256L * 2L)

## Data that does not fit is scaled across the whole of the type's range,
## which is what makes a narrow integer usable for continuous data at all
scaled <- asPacked(image, "int16")
expect_true(scaled@slope != 1)
expect_true(max(abs(as.array(scaled) - values)) < diff(range(values)) / 60000)
expect_equal(length(scaled@values), 256L * 2L)

## A quarter of the memory of a double image
expect_equal(length(scaled@values), length(values) * 8L / 4L)

## Every type round-trips within its own resolution
resolution <- c(int8 = 1 / 250, uint8 = 1 / 250, int16 = 1 / 60000,
                uint16 = 1 / 60000, int32 = 1e-8, float32 = 1e-6)
for (type in names(resolution))
{
    packed <- asPacked(image, type)
    expect_equal(storageType(packed), type, info = type)
    expect_equal(length(packed@values), 256L * unname(imply:::storageTypeSize[[type]]), info = type)
    expect_true(max(abs(as.array(packed) - values)) <= diff(range(values)) * resolution[[type]],
                info = paste(type, "round trip"))
}

## A constant image has no range to scale across, and must still round-trip
flat <- asPacked(denseImage(array(7, c(4L, 4L))), "int16")
expect_equal(as.array(flat), array(7, c(4L, 4L)))

expect_error(asPacked(image, "float64"), "one of")
expect_error(asPacked(denseImage(array(complex(real = 1, imaginary = 1), c(2L, 2L)))), "Only logical")

## --- Missing values --------------------------------------------------------

withNA <- values
withNA[5L] <- NA
naImage <- denseImage(withNA)

## Only a floating point type can carry missingness. Packing to an integer
## type would silently lose it, so it is refused instead
packedNA <- asPacked(naImage, "float32")
expect_true(is.na(as.array(packedNA)[5L]))
expect_equal(sum(is.na(as.array(packedNA))), 1L)
expect_error(asPacked(naImage, "int16"), "can only be packed as float32")
expect_error(asPacked(naImage, "uint8"), "can only be packed as float32")

## --- Summaries -------------------------------------------------------------

## Accumulation is in double whatever the storage type, so a summary of the
## packed image is exactly the summary of its unpacked values
unpacked <- as.array(f)
expect_identical(sum(f), sum(unpacked))
expect_identical(min(f), min(unpacked))
expect_identical(max(f), max(unpacked))
expect_identical(range(f), range(unpacked))
expect_identical(mean(f), mean(unpacked))

## ...and so agrees with the original within the type's precision
expect_equal(sum(f), sum(values), tolerance = 1e-6)
expect_equal(mean(f), mean(values), tolerance = 1e-6)

## Error does not compound over many elements, which it would if the
## accumulator were as narrow as the storage
constant <- asPacked(denseImage(array(0.1, c(100L, 100L, 100L))), "float32")
expect_equal(sum(constant), 1e5, tolerance = 1e-6)

## Missing values behave as they do elsewhere in R
expect_true(is.na(sum(packedNA)))
expect_equal(sum(packedNA, na.rm = TRUE), sum(withNA, na.rm = TRUE), tolerance = 1e-6)
expect_equal(mean(packedNA, na.rm = TRUE), mean(withNA, na.rm = TRUE), tolerance = 1e-6)

## --- Indexing --------------------------------------------------------------

## Reading a few values does not materialise the image
expect_equal(f[5L], values[5L], tolerance = 1e-6)
expect_equal(f[c(1L, 50L, 200L)], values[c(1L, 50L, 200L)], tolerance = 1e-6)
expect_equal(f[], as.array(f))
expect_equal(f[, , 1L], values[, , 1L], tolerance = 1e-6)

locations <- matrix(c(1L, 1L, 1L, 8L, 8L, 4L), ncol = 3L, byrow = TRUE)
expect_equal(f[locations], values[locations], tolerance = 1e-6)

## --- Applying functions ----------------------------------------------------

## The gather widens narrow values to double, so the function applied never
## learns that the image was stored narrowly
series <- denseImage(array(rnorm(6L * 7L * 8L * 5L), c(6L, 7L, 8L, 5L)))
packedSeries <- asPacked(series, "float32")
reference <- as.array(packedSeries)

expect_equal(imapply(packedSeries, 4, sum), apply(reference, 4, sum))
expect_equal(as.array(voxelApply(packedSeries, mean)), apply(reference, 1:3, mean))
expect_equal(lineApply(packedSeries, sum, axis = 1), apply(reference, c(2, 3), sum))
expect_equal(sliceApply(packedSeries, sum, axis = 3), apply(reference, 3, sum))

expect_equal(unique(as.vector(voxelApply(packedSeries, typeof))), "double")
expect_equal(unique(as.vector(voxelApply(packedSeries, length))), 5L)

## voxelApply on a packed image still yields an image, with the geometry
expect_true(isDenseImage(voxelApply(packedSeries, mean)))

## Dividing the work changes nothing
expect_identical(imapply(packedSeries, 4, sum, threads = 2L),
                 imapply(packedSeries, 4, sum, threads = 1L))

## --- Sparse images go through the same loop --------------------------------

## The generalisation that lets a packed image be applied over also lets a
## sparse one be, with an absent location becoming a zero during the gather
thresholded <- replace(as.array(series), abs(as.array(series)) < 1, 0)
sparseSeries <- asSparse(denseImage(thresholded))

expect_identical(imapply(sparseSeries, 4, sum), apply(thresholded, 4, sum))
expect_identical(as.array(voxelApply(sparseSeries, mean)), apply(thresholded, 1:3, mean))
expect_identical(lineApply(sparseSeries, sum, axis = 2), apply(thresholded, c(1, 3), sum))
expect_identical(imapply(sparseSeries, 4, sum, threads = 2L),
                 imapply(sparseSeries, 4, sum, threads = 1L))

## Nothing is materialised: the function sees the same values either way
expect_identical(imapply(sparseSeries, 4, range), apply(thresholded, 4, range))

## --- Arithmetic ------------------------------------------------------------

## Arithmetic on a packed image materialises, because keeping the result
## packed would need a scaling chosen from the result, which has to be
## computed first anyway
sum2 <- f + 1
expect_true(isDenseImage(sum2))
expect_equal(as.array(sum2), unpacked + 1)
expect_equal(voxelSize(sum2), c(2, 2, 2))

expect_true(isDenseImage(f * 2))
expect_equal(as.array(f * 2), unpacked * 2)
expect_true(isDenseImage(2 * f))
expect_equal(as.array(2 * f), 2 * unpacked)
expect_true(isDenseImage(f - f))
expect_equal(as.array(f - f), unpacked - unpacked)

## Mixing with a dense image works in both directions
expect_equal(as.array(f * image), unpacked * values)
expect_equal(as.array(image * f), values * unpacked)

## Repacking the answer is one call away
expect_true(isPackedImage(asPacked(f + 1, "float32")))

## --- Conversions between representations -----------------------------------

expect_identical(asPacked(f, "float32"), f)
expect_true(isPackedImage(asPacked(asSparse(image), "float32")))
expect_equal(as.array(asPacked(asSparse(image), "float32")), values, tolerance = 1e-6)
expect_true(isDenseImage(asDenseImage(as.array(f))))

## --- Storage claim ---------------------------------------------------------

## The headline: a float32 image is half the size of the double one, and an
## int16 image a quarter
big <- denseImage(array(rnorm(64L * 64L * 64L), c(64L, 64L, 64L)))
doubleBytes <- as.numeric(object.size(as.array(big)))

expect_true(as.numeric(object.size(asPacked(big, "float32")@values)) < doubleBytes * 0.55)
expect_true(as.numeric(object.size(asPacked(big, "int16")@values)) < doubleBytes * 0.3)
expect_true(as.numeric(object.size(asPacked(big, "uint8")@values)) < doubleBytes * 0.15)
