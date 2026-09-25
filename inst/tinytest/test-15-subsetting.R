## Indexing packed and sparse images, in both directions. The reference is
## always base R indexing of the equivalent dense array: whatever an
## expression does to the array, it must do to the image.

set.seed(15)
## Equal spatial extents, so that every subscript below is valid for every
## reorientation of the image too
values <- array(round(rnorm(5 * 5 * 5 * 2) * 3), c(5L, 5L, 5L, 2L))
values[abs(values) < 2] <- 0
las <- diag(c(-2, 2, 2.5, 1))
dense <- denseImage(values, worldTransform = las)

## Integer-valued data in int16 with no scaling, so packed values are exact
packed <- asPacked(dense, "int16")
expect_equal(packed@slope, 1)

## Images to test, each with the array it must behave like: plain, and
## reoriented so that indices pass through a non-trivial view
images <- list(
    sparse = list(image = asSparse(dense), reference = values),
    packed = list(image = packed, reference = values),
    sparseView = list(image = reorient(asSparse(dense), c(3L, -1L, 2L)), reference = as.array(reorient(dense, c(3L, -1L, 2L)))),
    packedView = list(image = reorient(packed, c(-2L, 3L, 1L)), reference = as.array(reorient(dense, c(-2L, 3L, 1L)))))

## --- Extraction ------------------------------------------------------------

extractions <- list(
    quote(X[2, 3, 4, 1]),
    quote(X[2, , 3, 1]),
    quote(X[, 2:3, -1, 2]),
    quote(X[c(TRUE, FALSE), , 5, ]),
    quote(X[1, 1, 1, 1, drop = FALSE]),
    quote(X[, , 1, 1, drop = FALSE]),
    quote(X[c(4, 1, 4), 5, , 2]),
    quote(X[]),
    quote(X[17]),
    quote(X[c(1, 99, 240, 99)]),
    quote(X[-(1:235)]),
    quote(X[cbind(c(1, 4), c(2, 5), c(5, 1), c(2, 1))]),
    quote(X[X > 2]),
    quote(X[integer(0)]),
    quote(X[NULL]))

for (name in names(images)) for (expression in extractions)
{
    info <- paste(name, deparse(expression))
    image <- images[[name]]$image
    reference <- images[[name]]$reference
    expected <- eval(expression, list(X = reference))
    result <- eval(expression, list(X = image))
    expect_equal(result, if (startsWith(name, "packed")) expected * 1 else expected, info = info)
    expect_identical(dim(result), dim(expected), info = info)
}

## --- Replacement -----------------------------------------------------------

replacements <- list(
    quote(X[2, 3, 4, 1] <- 7),
    quote(X[2, 3, 4, ] <- c(9, -9)),
    quote(X[, 2, , 1] <- 0),
    quote(X[, , , 2] <- 0),
    quote(X[1, , 1, ] <- 1:10),
    quote(X[c(TRUE, FALSE), -1, 2:3, 2] <- 5),
    quote(X[17] <- 11),
    quote(X[c(3, 3, 3)] <- c(1, 2, 3)),
    quote(X[X > 2] <- 0),
    quote(X[X == 0] <- 1),
    quote(X[cbind(c(1, 4), c(2, 5), c(5, 1), c(2, 1))] <- c(12, 13)),
    quote(X[] <- 0),
    quote(X[integer(0)] <- 5))

for (name in names(images)) for (expression in replacements)
{
    info <- paste(name, deparse(expression))
    image <- images[[name]]$image
    environment <- list2env(list(X = image))
    eval(expression, environment)
    result <- environment$X
    reference <- local({ X <- images[[name]]$reference; eval(expression); X })

    expect_identical(class(result), class(image), info = info)
    expect_equal(as.array(result), reference * 1, info = info)
    expect_identical(geometry(result), geometry(image), info = info)
    expect_identical(storageLayout(result), storageLayout(image), info = info)

    ## The original is untouched: replacement copies, as for any R object
    expect_equal(as.array(image), images[[name]]$reference * 1, info = info)
}

## A sparse image's mask always says where the data are: a location that
## gains a value is added, and one left with nothing but zeros is dropped
sparse <- asSparse(dense)
absent <- which(!mask(sparse))[1L]
location <- arrayInd(absent, dim(mask(sparse)))
sparse[location[1], location[2], location[3], 1] <- 4
expect_true(mask(sparse)[absent])
expect_equal(sum(mask(sparse)), sum(mask(asSparse(dense))) + 1)
sparse[location[1], location[2], location[3], 1] <- 0
expect_false(mask(sparse)[absent])
expect_identical(as.array(sparse), values)

present <- which(mask(sparse))[1L]
location <- arrayInd(present, dim(mask(sparse)))
sparse[location[1], location[2], location[3], ] <- 0
expect_false(mask(sparse)[present])
expect_equal(sum(mask(sparse)), sum(mask(asSparse(dense))) - 1)

## Missing values are not zero, so they are stored
sparse[1, 1, 1, 1] <- NA
expect_true(mask(sparse)[1, 1, 1])
expect_true(is.na(sparse[1, 1, 1, 1]))

## Values are promoted to a wider replacement type, as an array's are
counts <- asSparse(denseImage(array(c(0L, 1L, 0L, 2L), c(2L, 2L))))
expect_identical(typeof(counts@values), "integer")
counts[1, 1] <- 2.5
expect_identical(typeof(counts@values), "double")
expect_identical(as.array(counts), matrix(c(2.5, 1, 0, 2), 2L))
flags <- asSparse(denseImage(array(c(FALSE, TRUE), c(2L, 1L))))
flags[1, 1] <- 3L
expect_identical(as.array(flags), matrix(c(3L, 1L), 2L))

## --- Packed values and their scaling ---------------------------------------

## Values are rounded to the type's resolution under the existing scaling
scaled <- asPacked(denseImage(array(seq(0, 1, length.out = 24), c(2L, 3L, 4L))), "uint8")
scaled[1] <- 0.5
expect_equal(scaled[1], 0.5, tolerance = scaled@slope)
expect_equal(scaled@slope, asPacked(denseImage(array(seq(0, 1, length.out = 24), c(2L, 3L, 4L))), "uint8")@slope)

## ...and values it cannot reach, or missing values in an integer type, are
## refused rather than clamped or rescaled
expect_error(scaled[1] <- 2, "current scaling")
expect_error(packed[1] <- 1e6, "current scaling")
expect_error(packed[1] <- NA, "Missing values cannot be stored")
expect_error(packed[1] <- 1i, "Complex values")
floating <- asPacked(dense, "float32")
floating[2] <- NA
expect_true(is.na(floating[2]))

## An int32 never stores R's integer NA
wide <- asPacked(dense, "int32")
expect_error(wide[1] <- -2147483648, "current scaling")
wide[1] <- -2147483647
expect_equal(wide[1], -2147483647)

## --- Subscripts ------------------------------------------------------------

for (image in list(asSparse(dense), packed))
{
    expect_error(image[6, 1, 1, 1], "out of bounds")
    expect_error(image[251], "out of bounds")
    expect_error(image[c(1, -1)], "cannot be mixed")
    expect_error(image[NA_real_], "Missing values")
    expect_error(image[1, 1], "4 dimensions, but 2 subscripts")
    expect_error(image["a"], "by number or by logical value")
    expect_error(image[rep(TRUE, 251)], "too long")
    expect_error(image[1, 1, 1, 1] <- numeric(0), "length zero")
    expect_warning(image[1:3] <- 1:2, "not a multiple")
}

## Subscripts are evaluated where the call is made, as base R's are
setCorner <- function (image, value)
{
    i <- 1
    j <- 2:3
    image[i, j, 1, 1] <- value
    image[i, j, 1, 1]
}
expect_equal(setCorner(asSparse(dense), 8), c(8, 8))
expect_equal(setCorner(packed, -8), c(-8, -8))
