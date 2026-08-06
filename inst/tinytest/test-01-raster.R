## Raster geometry, indexing, line traversal, blocking and permuted views.
##
## The fixed-dimensionality and runtime-dimensionality code paths are the same
## algorithm instantiated two ways, so every test that can be run on one is run
## on both and the results compared.

## These are internal primitives for now. The user-facing verbs (imapply,
## voxelApply, sliceApply) will be layered on them, and exported in their place
rasterInfo <- imply:::rasterInfo
flattenIndices <- imply:::flattenIndices
expandIndices <- imply:::expandIndices
lineSums <- imply:::lineSums
blockPartition <- imply:::blockPartition
permuteView <- imply:::permuteView

dims <- c(4L, 5L, 3L, 2L)
x <- array(as.double(seq_len(prod(dims))), dims)

## --- Raster geometry -------------------------------------------------------

info <- rasterInfo(x)
expect_equal(info$dim, dims)
expect_equal(info$nDims, 4L)
expect_equal(info$size, prod(dims))
expect_true(info$contiguous)
expect_true(info$fixed)

## Strides are cumulative products with the first dimension moving fastest.
## They come back as double rather than integer deliberately: a stride can
## exceed integer range on a long vector
expect_equal(info$strides, cumprod(c(1, dims[-length(dims)])))
expect_true(is.double(info$strides))

## Leading three dimensions are spatial by default, so the trailing one indexes
## the value held at each location
expect_equal(info$spatial, 3L)
expect_equal(info$spatialSize, prod(dims[1:3]))
expect_equal(info$elementSize, dims[4])

## An explicit split is honoured, including the degenerate ends
expect_equal(rasterInfo(x, spatial = 2L)$elementSize, prod(dims[3:4]))
expect_equal(rasterInfo(x, spatial = 0L)$spatialSize, 1)
expect_equal(rasterInfo(x, spatial = 0L)$elementSize, prod(dims))
expect_equal(rasterInfo(x, spatial = 4L)$elementSize, 1)
expect_error(rasterInfo(x, spatial = 5L), "exceeds the dimensionality")

## The fixed and dynamic paths must agree in every respect but which they are
fixedInfo <- rasterInfo(x)
dynamicInfo <- rasterInfo(x, forceDynamic = TRUE)
shared <- setdiff(names(fixedInfo), "fixed")
expect_identical(fixedInfo[shared], dynamicInfo[shared])
expect_false(dynamicInfo$fixed)

## A plain vector with no dim attribute is one-dimensional, not an error
v <- as.double(1:10)
expect_equal(rasterInfo(v)$nDims, 1L)
expect_equal(rasterInfo(v)$dim, 10L)
expect_equal(rasterInfo(v)$spatial, 1L)

## Beyond five dimensions there is no fixed instantiation, so it falls through
big <- array(1, c(2L, 2L, 2L, 2L, 2L, 2L))
expect_false(rasterInfo(big)$fixed)
expect_equal(rasterInfo(big)$nDims, 6L)

## Five dimensions is the last fixed case
expect_true(rasterInfo(array(1, rep(2L, 5)))$fixed)

## --- Index flattening ------------------------------------------------------

## An array of consecutive integers indexed by a location matrix yields exactly
## the linear index of each location, which is the reference to check against
reference <- array(seq_len(prod(dims)), dims)
locs <- as.matrix(expand.grid(1:4, 1:5, 1:3, 1:2))
dimnames(locs) <- NULL

expect_equal(as.integer(flattenIndices(x, locs)), as.integer(reference[locs]))
expect_equal(flattenIndices(x, locs), flattenIndices(x, locs, forceDynamic = TRUE))

## And the inverse
expect_equal(expandIndices(x, seq_len(prod(dims))), arrayInd(seq_len(prod(dims)), dims))
expect_equal(expandIndices(x, seq_len(prod(dims))),
             expandIndices(x, seq_len(prod(dims)), forceDynamic = TRUE))

## Round trip
expect_equal(as.integer(flattenIndices(x, expandIndices(x, seq_len(prod(dims))))),
             seq_len(prod(dims)))

expect_error(flattenIndices(x, matrix(c(5L, 1L, 1L, 1L), nrow = 1)), "out of range")
expect_error(flattenIndices(x, matrix(1L, nrow = 1, ncol = 3)), "3 columns")

## --- Line traversal --------------------------------------------------------

## Summing along every line in a direction is the complement of applying sum
## over the other margins, and lines are enumerated in the same order
for (d in 1:4)
    expect_equal(lineSums(x, d), as.vector(apply(x, (1:4)[-d], sum)),
                 info = paste("lineSums along dimension", d))

for (d in 1:4)
    expect_equal(lineSums(x, d), lineSums(x, d, forceDynamic = TRUE))

## Integer and logical input are accepted and accumulate in double
xi <- array(seq_len(prod(dims)), dims)
expect_equal(lineSums(xi, 2), as.vector(apply(xi, (1:4)[-2], sum)))
xl <- array(rep(c(TRUE, FALSE), length.out = prod(dims)), dims)
expect_equal(lineSums(xl, 1), as.vector(apply(xl, (1:4)[-1], sum)))

## Missing values propagate the way sum() does, for both real and integer input
xna <- x
xna[2, 2, 2, 2] <- NA
expect_equal(lineSums(xna, 1), as.vector(apply(xna, (1:4)[-1], sum)))
expect_true(anyNA(lineSums(xna, 1)))

xina <- xi
xina[2, 2, 2, 2] <- NA
expect_equal(lineSums(xina, 1), as.vector(apply(xina, (1:4)[-1], sum)))

expect_error(lineSums(x, 5), "out of range")
expect_error(lineSums(array(complex(real = 1, imaginary = 1), c(2, 2)), 1), "not yet supported")

## --- Blocking --------------------------------------------------------------

## Blocks must tile the spatial locations exactly once, in order
checkPartition <- function (partition, total) {
    starts <- unname(partition[, "start"])
    lengths <- unname(partition[, "length"])
    expect_equal(starts[1], 1)
    expect_equal(sum(lengths), total)
    expect_equal(starts[-1], head(starts + lengths, -1))
}

checkPartition(blockPartition(x, targetElements = 8), prod(dims[1:3]))
checkPartition(blockPartition(x, targetElements = 1), prod(dims[1:3]))
checkPartition(blockPartition(x, targetElements = 1e6), prod(dims[1:3]))
expect_equal(nrow(blockPartition(x, targetElements = 1e6)), 1L)

## Splitting into a fixed number of pieces is what makes a requested thread
## count meaningful, so it must never exceed the number asked for
for (n in c(1L, 2L, 3L, 8L, 100L)) {
    partition <- blockPartition(x, count = n)
    checkPartition(partition, prod(dims[1:3]))
    expect_true(nrow(partition) <= max(n, 1L), info = paste("count =", n))
}

## --- Permuted views --------------------------------------------------------

## Permutation reorders the stride vector rather than moving data, so reading
## through it must reproduce aperm() exactly
for (perm in list(c(1, 2, 3, 4), c(4, 3, 2, 1), c(2, 1, 4, 3), c(3, 1, 4, 2))) {
    expect_equal(permuteView(x, perm), aperm(x, perm),
                 info = paste("permutation", paste(perm, collapse = "")))
    expect_equal(permuteView(x, perm), permuteView(x, perm, forceDynamic = TRUE))
}

## Including for the dynamic path beyond five dimensions
expect_equal(permuteView(big, c(6, 5, 4, 3, 2, 1)), aperm(big, c(6, 5, 4, 3, 2, 1)))

## And for every storage mode
expect_equal(permuteView(xi, c(2, 1, 3, 4)), aperm(xi, c(2, 1, 3, 4)))
expect_equal(permuteView(xl, c(2, 1, 3, 4)), aperm(xl, c(2, 1, 3, 4)))
xc <- array(complex(real = seq_len(prod(dims)), imaginary = 1), dims)
expect_equal(permuteView(xc, c(2, 1, 3, 4)), aperm(xc, c(2, 1, 3, 4)))

expect_error(permuteView(x, c(1, 2, 3)), "length 3")
expect_error(permuteView(x, c(1, 2, 3, 5)), "out of range")
expect_error(permuteView(x, c(1, 1, 2, 3)), "valid ordering")

## --- Zero-copy -------------------------------------------------------------

## Read-only entry points must borrow R's memory rather than duplicating it,
## and must not coerce the input in passing
address <- imply:::dataAddress(x)
invisible(lineSums(x, 1))
invisible(rasterInfo(x))
invisible(permuteView(x, c(2, 1, 3, 4)))
expect_identical(imply:::dataAddress(x), address)

## An integer array must not be silently promoted to double on the way in
addressInt <- imply:::dataAddress(xi)
invisible(lineSums(xi, 1))
expect_identical(imply:::dataAddress(xi), addressInt)
expect_true(is.integer(xi))

## --- Degenerate shapes -----------------------------------------------------

empty <- array(numeric(0), c(0L, 3L, 2L))
expect_equal(rasterInfo(empty)$size, 0)
expect_equal(nrow(blockPartition(empty)), 0L)
expect_equal(length(permuteView(empty, c(3, 2, 1))), 0L)

unitary <- array(1, c(1L, 1L, 1L))
expect_equal(permuteView(unitary, c(2, 3, 1)), aperm(unitary, c(2, 3, 1)))
expect_equal(lineSums(unitary, 1), 1)
