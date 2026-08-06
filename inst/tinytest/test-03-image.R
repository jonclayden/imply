## The denseImage class: construction, validation, geometry and coercion.

## --- Construction ----------------------------------------------------------

x <- denseImage(array(as.double(1:24), c(2L, 3L, 4L)))
expect_true(isDenseImage(x))
expect_true(is.array(x))
expect_equal(dim(x), c(2L, 3L, 4L))
expect_equal(spatial(x), 3L)
expect_equal(pixdim(x), c(1, 1, 1))
expect_equal(xform(x), diag(4))
expect_equal(orientation(x), "RAS")

## Every storage mode that occurs in practice must survive, masks being
## logical and phase data complex. This is why the class does not inherit from
## an S7 base type, which would pin it to one mode
for (mode in c("logical", "integer", "double", "complex")) {
    a <- array(as.vector(1:24, mode = mode), c(2L, 3L, 4L))
    image <- denseImage(a)
    expect_identical(typeof(image), mode, info = mode)
    expect_true(isDenseImage(image), info = mode)
    expect_equal(dim(image), c(2L, 3L, 4L), info = mode)
}

expect_error(denseImage(array("a", c(2, 2))), "logical, integer, double or complex")

## A bare vector becomes one-dimensional rather than an error
v <- denseImage(1:10)
expect_equal(dim(v), 10L)
expect_equal(spatial(v), 1L)
expect_equal(pixdim(v), 1)

## --- The spatial/element split ---------------------------------------------

## A four-dimensional image is three spatial dimensions with a time series at
## each location by default
series <- denseImage(array(0, c(4L, 4L, 4L, 10L)))
expect_equal(spatial(series), 3L)
expect_equal(length(pixdim(series)), 3L)

## ...but the split can be stated explicitly, for a vector field say
field <- denseImage(array(0, c(4L, 4L, 4L, 3L)), spatial = 3L, pixdim = c(2, 2, 2))
expect_equal(spatial(field), 3L)
slices <- denseImage(array(0, c(4L, 4L, 6L)), spatial = 2L, pixdim = c(1, 1))
expect_equal(spatial(slices), 2L)
expect_equal(pixdim(slices), c(1, 1))

expect_error(denseImage(array(0, c(2, 2)), spatial = 3L), "between 0 and 2")
expect_error(denseImage(array(0, c(2, 2, 2)), pixdim = c(1, 1)), "one element per spatial dimension")

## The compiled side reads the split straight off the attribute, so the class
## and the raster agree without any S7 knowledge in C++
expect_equal(imply:::rasterInfo(slices)$spatial, 2L)
expect_equal(imply:::rasterInfo(slices)$elementSize, 6)
expect_equal(imply:::rasterInfo(series)$elementSize, 10)
expect_equal(imply:::rasterInfo(x)$spatial, 3L)

## --- Validation ------------------------------------------------------------

expect_error(denseImage(array(0, c(2, 2, 2)), xform = diag(3)), "4x4")
expect_error(denseImage(array(0, c(2, 2, 2)), xform = rbind(diag(4)[1:3, ], c(1, 1, 1, 1))), "affine")
expect_error(denseImage(array(0, c(2, 2, 2)), pixdim = c(1, NA, 1)), "must not be missing")

## --- Geometry accessors ----------------------------------------------------

las <- rbind(c(-2, 0, 0, 90), c(0, 2, 0, -126), c(0, 0, 2, -72), c(0, 0, 0, 1))
image <- denseImage(array(0, c(10L, 10L, 10L)), pixdim = c(2, 2, 2), xform = las)
expect_equal(xform(image), las)
expect_equal(orientation(image), "LAS")

## Replacing the voxel dimensions must carry the transform with them, or the
## two would disagree
pixdim(image) <- c(3, 3, 3)
expect_equal(pixdim(image), c(3, 3, 3))
expect_equal(diag(xform(image))[1:3], c(3, 3, 3))
expect_true(isDenseImage(image))

xform(image) <- las
expect_equal(xform(image), las)
expect_error(xform(image) <- diag(3), "4x4")
expect_error(pixdim(image) <- c(1, 1), "one element per spatial dimension")

## A template supplies whatever is not given explicitly
copy <- denseImage(array(1, c(10L, 10L, 10L)), template = image)
expect_equal(xform(copy), xform(image))
expect_equal(pixdim(copy), pixdim(image))

## --- S7 properties ---------------------------------------------------------

## Properties are readable through @, and are stored as ordinary attributes so
## that compiled code needs no knowledge of S7
expect_equal(x@spatial, 3L)
expect_equal(x@pixdim, c(1, 1, 1))
expect_equal(attr(x, "pixdim"), x@pixdim)
expect_equal(attr(x, "spatial"), x@spatial)

## S7 qualifies the class name with the package, so a bare inherits() test
## fails. This is exactly why isDenseImage() exists
expect_equal(class(x)[1], "imply::denseImage")
expect_false(inherits(x, "denseImage"))
expect_true(isDenseImage(x))

## The validator runs on assignment, not merely at construction
bad <- x
expect_error({bad@pixdim <- c(1, 1)}, "one element per spatial dimension")
expect_error({bad@xform <- diag(3)}, "4x4")
expect_error({bad@spatial <- 9L}, "between 0 and 3")
expect_error({bad@pixdim <- c(1, NA, 1)}, "must not be missing")

## ...and a valid assignment goes through
good <- x
good@pixdim <- c(2, 2, 2)
expect_equal(good@pixdim, c(2, 2, 2))
expect_true(isDenseImage(good))

## --- Subsetting ------------------------------------------------------------

## S7 objects are not subsettable by default, so these methods are what make
## an image indexable at all. The result is a plain array: an arbitrary index
## has no well-defined geometry
expect_equal(x[1, 1, 1], 1)
expect_equal(dim(x[, , 1]), c(2L, 3L))
expect_false(isDenseImage(x[, , 1]))
expect_equal(x[, , 1], as.array(x)[, , 1])
expect_equal(x[x > 20], c(21, 22, 23, 24))
expect_equal(dim(x[, , 1, drop = FALSE]), c(2L, 3L, 1L))
expect_equal(as.vector(x[]), as.vector(as.array(x)))

## Replacement preserves shape, so the geometry still applies and the result
## is still an image
replaced <- x
replaced[1, 1, 1] <- 99
expect_true(isDenseImage(replaced))
expect_equal(replaced[1, 1, 1], 99)
expect_equal(pixdim(replaced), pixdim(x))
expect_equal(xform(replaced), xform(x))
expect_equal(dim(replaced), dim(x))

## --- Coercion --------------------------------------------------------------

plain <- as.array(x)
expect_false(isDenseImage(plain))
expect_null(attr(plain, "pixdim"))
expect_null(attr(plain, "xform"))
expect_null(attr(plain, "spatial"))
expect_equal(dim(plain), dim(x))
expect_equal(as.vector(plain), as.vector(x))

## asDenseImage leaves an image alone but promotes a plain array
expect_identical(asDenseImage(x), x)
expect_true(isDenseImage(asDenseImage(array(0, c(2, 2)))))

## --- Base R behaviour is inherited -----------------------------------------

## Arithmetic and summaries work with no methods of their own, and arithmetic
## keeps the geometry
expect_equal(sum(x), sum(1:24))
expect_equal(max(x), 24)
expect_true(isDenseImage(x * 2))
expect_equal(pixdim(x * 2), pixdim(x))
expect_equal(as.vector(x * 2), as.vector(1:24) * 2)
expect_true(isDenseImage(sqrt(x)))

## Plain arrays remain acceptable everywhere: no custom class is ever required
expect_equal(imply:::lineSums(as.array(x), 1), imply:::lineSums(x, 1))
expect_equal(orientation(diag(4)), orientation(denseImage(array(0, c(2, 2, 2)))))

## --- Printing --------------------------------------------------------------

output <- capture.output(print(x))
expect_true(any(grepl("Dense image: 2 x 3 x 4", output)))
expect_true(any(grepl("RAS", output)))
expect_true(any(grepl("double", output)))
expect_true(any(grepl("Values per location: 10", capture.output(print(series)))))
