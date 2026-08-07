## Mask-aware application, and the packed values as a matrix.

set.seed(1)
spatialDims <- c(8L, 7L, 6L)
nT <- 12L
brain <- array(runif(prod(spatialDims)) < 0.4, spatialDims)
values <- array(rnorm(prod(spatialDims) * nT), c(spatialDims, nT))
values[!brain] <- 0

image <- denseImage(values, pixdim = c(2, 2, 2))
packed <- asSparse(image)
selected <- sum(brain)

## --- maskedMatrix ----------------------------------------------------------

## Packing is voxel-major, so the stored values already are the masked data
## matrix: one column per stored location, values at a location contiguous
m <- maskedMatrix(packed)
expect_equal(dim(m), c(nT, selected))
expect_true(is.matrix(m))

## Each column is the series at the corresponding stored location
stored <- which(as.vector(mask(packed)))
for (k in c(1L, 2L, length(stored)))
{
    location <- arrayInd(stored[k], spatialDims)
    expect_equal(as.vector(m[, k]),
                 as.vector(values[location[1L], location[2L], location[3L], ]),
                 info = paste("column", k))
}

## It costs nothing, because the values are held in this shape already
expect_identical(imply:::dataAddress(m), imply:::dataAddress(packed@values))

## A scalar image gives one row
scalar <- asSparse(denseImage(array(c(1, 0, 0, 2, 0, 3), c(3L, 2L)), spatial = 2L, pixdim = c(1, 1)))
expect_equal(dim(maskedMatrix(scalar)), c(1L, 3L))
expect_equal(as.vector(maskedMatrix(scalar)), c(1, 2, 3))

## An empty image has no columns
empty <- asSparse(denseImage(array(0, c(4L, 4L, 4L))))
expect_equal(ncol(maskedMatrix(empty)), 0L)

expect_error(maskedMatrix(image), "Only a sparse image")

## --- Mask-aware voxelApply -------------------------------------------------

## The point of the mask: locations outside it are not visited at all
calls <- 0L
invisible(voxelApply(image, function (v) { calls <<- calls + 1L; 0 }, mask = brain))
expect_equal(calls, selected)
expect_true(calls < prod(spatialDims))

## Without a mask every location is visited, which is what makes the mask worth
## having rather than relying on the image being sparse
calls <- 0L
invisible(voxelApply(packed, function (v) { calls <<- calls + 1L; 0 }))
expect_equal(calls, prod(spatialDims))

## Answers inside the mask match the unmasked ones, and outside is fill
reference <- voxelApply(image, function (v) mean(v))
masked <- voxelApply(image, function (v) mean(v), mask = brain)

expect_true(isDenseImage(masked))
expect_equal(dim(masked), spatialDims)
expect_equal(as.array(masked)[brain], as.array(reference)[brain])
expect_true(all(as.array(masked)[!brain] == 0))

## The fill value is settable, including to NA
withNA <- voxelApply(image, function (v) mean(v), mask = brain, fill = NA)
expect_true(all(is.na(as.array(withNA)[!brain])))
expect_equal(as.array(withNA)[brain], as.array(reference)[brain])
expect_equal(as.array(voxelApply(image, function (v) mean(v), mask = brain, fill = -1))[!brain][1L], -1)

## Geometry is carried through
expect_equal(pixdim(masked), c(2, 2, 2))
expect_equal(xform(masked), xform(image))

## --- Results wider than one value per location -----------------------------

## A vector result gains a leading dimension, as it does unmasked, and the
## fill covers every element at an absent location
ranges <- voxelApply(image, function (v) range(v), mask = brain)
expect_equal(dim(ranges), c(2L, spatialDims))
expect_false(isDenseImage(ranges))

unmaskedRanges <- voxelApply(image, function (v) range(v))
inside <- rep(as.vector(brain), each = 2L)
expect_equal(as.vector(ranges)[inside], as.vector(unmaskedRanges)[inside])
expect_true(all(as.vector(ranges)[!inside] == 0))

## A result the same length as the input round-trips to the original values
identity <- voxelApply(image, function (v) v, mask = brain)
expect_equal(dim(identity), c(nT, spatialDims))
expect_equal(as.vector(identity)[rep(as.vector(brain), each = nT)],
             as.vector(aperm(values, c(4, 1, 2, 3)))[rep(as.vector(brain), each = nT)])

## --- Where the mask comes from ---------------------------------------------

## A sparse image's own mask, given directly
expect_equal(voxelApply(image, function (v) mean(v), mask = packed),
             voxelApply(image, function (v) mean(v), mask = brain))

## A numeric mask, where non-zero selects
expect_equal(voxelApply(image, function (v) mean(v), mask = brain * 1.0),
             voxelApply(image, function (v) mean(v), mask = brain))

## Applying a sparse image's own mask to itself uses the stored values with no
## copy, and must give the same answer as the dense route
expect_equal(as.array(voxelApply(packed, function (v) mean(v), mask = brain)),
             as.array(voxelApply(image, function (v) mean(v), mask = brain)))

## A mask that is not the image's own still works, by gathering
half <- brain & (slice.index(brain, 3) <= 3)
expect_equal(sum(voxelApply(packed, function (v) 1, mask = half)), sum(half))

expect_error(voxelApply(image, mean, mask = array(TRUE, c(2L, 2L))), "one element per spatial location")
expect_error(voxelApply(image, mean, mask = array(NA, spatialDims)), "must not contain missing")
expect_error(voxelApply(image, mean, mask = array(FALSE, spatialDims)), "selects no locations")
expect_error(voxelApply(image, mean, mask = "brain"), "logical array")

## --- Interaction with the other arguments ----------------------------------

expect_identical(voxelApply(image, function (v) mean(v), mask = brain, threads = 2L),
                 voxelApply(image, function (v) mean(v), mask = brain, threads = 1L))

reported <- 0L
invisible(voxelApply(image, function (v) mean(v), mask = brain,
                     progress = function (done, total) reported <<- reported + 1L))
expect_true(reported > 0L)

## Extra arguments reach the function
withMissing <- values
withMissing[which(brain)[1L]] <- NA
expect_false(is.na(as.array(voxelApply(denseImage(withMissing), function (v) mean(v, na.rm = TRUE),
                                       mask = brain))[which(brain)[1L]]))

## A narrow image works too, its values widened on the way in
expect_equal(as.array(voxelApply(asPacked(image, "float32"), function (v) mean(v), mask = brain)),
             as.array(voxelApply(image, function (v) mean(v), mask = brain)),
             tolerance = 1e-6)
