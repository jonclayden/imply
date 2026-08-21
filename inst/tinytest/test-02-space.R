## Image geometry: affine inversion, voxel size/transform decomposition, and
## coordinate conversion.

## A typical clinical transform: 2 mm isotropic, left-handed x axis
las <- rbind(c(-2, 0, 0,   90),
             c( 0, 2, 0, -126),
             c( 0, 0, 2,  -72),
             c( 0, 0, 0,    1))

## --- Affine inversion ------------------------------------------------------

expect_equal(imply:::invertXform(las), solve(las))
expect_equal(imply:::invertXform(las) %*% las, diag(4))
expect_equal(imply:::invertXform(imply:::invertXform(las)), las)

## An oblique transform, to check the general path rather than the diagonal one
set.seed(1)
oblique <- diag(4)
oblique[1:3, 1:3] <- qr.Q(qr(matrix(rnorm(9), 3))) %*% diag(c(1.2, 0.8, 3))
oblique[1:3, 4] <- c(10, -20, 30)
expect_equal(imply:::invertXform(oblique), solve(oblique))

expect_error(imply:::invertXform(diag(3)), "4x4")
expect_error(imply:::invertXform(rbind(diag(4)[1:3, ], c(1, 1, 1, 1))), "affine")
expect_error(imply:::invertXform(diag(c(1, 1, 0, 1))), "singular")

## --- Decomposition into voxel size and world transform ----------------------

## las decomposes to voxel size (2, 2, 2) and a rigid, reflected frame
decomposedLas <- imply:::decomposeTransform(las, 3L)
expect_equal(decomposedLas$voxelSize, c(2, 2, 2))
expect_equal(decomposedLas$orientation %*% diag(c(2, 2, 2, 1)), las)

## Recomposing gets back the original affine
expect_equal(imply:::composeTransform(decomposedLas$orientation, decomposedLas$voxelSize), las)

## Anisotropic, oblique voxel size and orientation round-trip too
decomposedOblique <- imply:::decomposeTransform(oblique, 3L)
expect_equal(decomposedOblique$voxelSize, c(1.2, 0.8, 3))
expect_equal(imply:::composeTransform(decomposedOblique$orientation, decomposedOblique$voxelSize), oblique)

## A genuinely sheared affine -- as might arise from a 12-parameter affine
## registration to a template space -- cannot be decomposed into rotation and
## voxel size, and is rejected rather than silently mangled
sheared <- diag(4)
sheared[1:3, 1:3] <- matrix(c(1, 0, 0, 0.3, 1, 0, 0, 0, 1), 3)
expect_error(imply:::decomposeTransform(sheared, 3L), "shear")

## --- Voxel size and world transform never desync ----------------------------

image <- denseImage(array(0, c(91L, 109L, 91L)), voxelSize = c(2, 2, 2), worldTransform = las)
expect_equal(worldTransform(image), las)
expect_equal(voxelSize(image), c(2, 2, 2))

## Setting voxel size alone must leave rotation and translation untouched --
## this is the regression case for the original bug, where pixdim<- discarded
## both
voxelSize(image) <- c(3, 3, 3)
expect_equal(worldTransform(image)[1:3, 4], las[1:3, 4])
expect_equal(worldTransform(image)[1:3, 1:3], las[1:3, 1:3] / 2 * 3)
expect_equal(voxelSize(image), c(3, 3, 3))

## worldTransform<- decomposes its argument the same way and, symmetrically,
## must leave nothing of the previous voxel size behind
worldTransform(image) <- oblique
expect_equal(worldTransform(image), oblique)
expect_equal(voxelSize(image), c(1.2, 0.8, 3))

expect_error(worldTransform(image) <- sheared, "shear")

## --- Coordinate conversion -------------------------------------------------

image <- denseImage(array(0, c(91L, 109L, 91L)), voxelSize = c(2, 2, 2), worldTransform = las)

## Converting a world point to voxel coordinates must invert the transform.
## The implementation this was ported from applied the forward transform here,
## so a round trip is the regression test for that bug
voxel <- c(10, 20, 30)
world <- fromVoxel(voxel, image)
expect_equal(as.vector(world), as.vector(las %*% c(voxel, 1))[1:3])
expect_equal(as.vector(toVoxel(world, image)), voxel)

## Explicitly: the reverse conversion is not the forward transform
expect_false(isTRUE(all.equal(as.vector(toVoxel(world, image)),
                              as.vector(las %*% c(world, 1))[1:3])))

## Round trip over many points, including the oblique case
set.seed(7)
points <- matrix(runif(300, 0, 80), ncol = 3)
expect_equal(toVoxel(fromVoxel(points, image), image), points)

obliqueImage <- denseImage(array(0, c(10L, 10L, 10L)), voxelSize = c(1.2, 0.8, 3), worldTransform = oblique)
expect_equal(toVoxel(fromVoxel(points, obliqueImage), obliqueImage), points)

## Scaled coordinates apply the voxel dimensions but ignore rotation
expect_equal(as.vector(fromVoxel(c(1, 1, 1), image, type = "scaled")), c(2, 2, 2))
expect_equal(as.vector(toVoxel(c(2, 2, 2), image, type = "scaled")), c(1, 1, 1))

## Voxel coordinates pass through untouched
expect_equal(as.vector(fromVoxel(c(3, 4, 5), image, type = "voxel")), c(3, 4, 5))

expect_error(fromVoxel(matrix(1:8, ncol = 4), image), "three columns")

## --- Rounding --------------------------------------------------------------

points <- matrix(c(1.2, 2.7, 3.5, -0.4, 9.9, 4.5), ncol = 3, byrow = TRUE)

expect_equal(imply:::roundPoints(points, "none"), points)

## Conventional rounding means exactly what R's round() means, including
## breaking an exact half to even rather than away from zero
expect_equal(imply:::roundPoints(points, "conventional"), round(points))
expect_equal(as.vector(imply:::roundPoints(matrix(c(3.5, 4.5, -0.5), ncol = 3), "conventional")),
             c(4, 4, 0))

## Probabilistic rounding always lands on one of the two neighbouring integers
set.seed(1)
for (i in 1:20) {
    rounded <- imply:::roundPoints(points, "probabilistic")
    expect_true(all(rounded >= floor(points) & rounded <= ceiling(points)))
}

## ...and never steps beyond a stated bound
set.seed(2)
edge <- matrix(c(9.6, 9.6, 9.6), ncol = 3)
for (i in 1:50)
    expect_true(all(imply:::roundPoints(edge, "probabilistic", bounds = c(10, 10, 10)) <= 9))

## Randomness comes from the package's own Mersenne twister rather than R's
## global generator, so that it is safe to call from a worker thread. The seed
## is still drawn from R, so set.seed() must remain in charge of reproducibility
set.seed(11)
first <- imply:::roundPoints(points, "probabilistic")
set.seed(11)
expect_identical(imply:::roundPoints(points, "probabilistic"), first)

set.seed(12)
expect_false(identical(imply:::roundPoints(points, "probabilistic"), first))

## Over many draws it is genuinely stochastic, and biased toward the nearer
## integer rather than uniform
set.seed(3)
draws <- replicate(2000, imply:::roundPoints(matrix(c(5.25, 5.25, 5.25), ncol = 3), "probabilistic")[1])
expect_true(all(draws %in% c(5, 6)))
expect_true(mean(draws == 5) > 0.6 && mean(draws == 5) < 0.9)

expect_error(imply:::roundPoints(points, "nonsense"), "Rounding type")
expect_error(imply:::pointsToVoxel(points, diag(4), c(1, 1, 1), "nonsense"), "Point type")
