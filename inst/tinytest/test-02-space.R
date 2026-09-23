## Image geometry: voxel size/transform decomposition and coordinate
## conversion. Affine inversion has no dedicated export or test of its own --
## Affine::inverse() is exercised indirectly, but thoroughly, by the
## world-type round trips below, since it is only ever handed an orthonormal
## (and so always invertible) matrix by any path reachable from R

## A typical clinical transform: 2 mm isotropic, left-handed x axis
las <- rbind(c(-2, 0, 0,   90),
             c( 0, 2, 0, -126),
             c( 0, 0, 2,  -72),
             c( 0, 0, 0,    1))

## An oblique transform, to check the general path rather than the diagonal one
set.seed(1)
oblique <- diag(4)
oblique[1:3, 1:3] <- qr.Q(qr(matrix(rnorm(9), 3))) %*% diag(c(1.2, 0.8, 3))
oblique[1:3, 4] <- c(10, -20, 30)

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

## The setters must preserve whatever image class they are given rather than
## densifying as a side effect -- this is the regression case for a bug where
## they reached for asDenseImage(), which does not understand sparse or
## packed images and so refused them outright
sparseImg <- asSparse(image)
voxelSize(sparseImg) <- c(4, 4, 4)
expect_true(isSparseImage(sparseImg))
expect_equal(voxelSize(sparseImg), c(4, 4, 4))
expect_equal(worldTransform(sparseImg)[1:3, 4], oblique[1:3, 4])

worldTransform(sparseImg) <- las
expect_true(isSparseImage(sparseImg))
expect_equal(worldTransform(sparseImg), las)

packedImg <- asPacked(image, "int16")
voxelSize(packedImg) <- c(5, 5, 5)
expect_true(isPackedImage(packedImg))
expect_equal(voxelSize(packedImg), c(5, 5, 5))

worldTransform(packedImg) <- las
expect_true(isPackedImage(packedImg))
expect_equal(worldTransform(packedImg), las)

## --- Coordinate conversion -------------------------------------------------

image <- denseImage(array(0, c(91L, 109L, 91L)), voxelSize = c(2, 2, 2), worldTransform = las)

## toVoxel()/fromVoxel() use one-based voxel coordinates, matching x[i,j,k],
## even though the affine itself is zero-based throughout. The first voxel
## therefore sits exactly at the affine's own translation, with no offset
## left for the caller to apply by hand
expect_equal(as.vector(fromVoxel(c(1, 1, 1), image)), las[1:3, 4])
expect_equal(as.vector(toVoxel(las[1:3, 4], image)), c(1, 1, 1))

## Converting a world point to voxel coordinates must invert the transform.
## The implementation this was ported from applied the forward transform here,
## so a round trip is the regression test for that bug. The affine is applied
## to the zero-based equivalent of the one-based input
voxel <- c(10, 20, 30)
world <- fromVoxel(voxel, image)
expect_equal(as.vector(world), as.vector(las %*% c(voxel - 1, 1))[1:3])
expect_equal(as.vector(toVoxel(world, image)), voxel)

## Explicitly: the reverse conversion is not the forward transform
expect_false(isTRUE(all.equal(as.vector(toVoxel(world, image)),
                              as.vector(las %*% c(world, 1))[1:3])))

## Round trip over many points, including the oblique case. This holds
## regardless of the one-based shift, since it cancels out in a round trip
set.seed(7)
points <- matrix(runif(300, 0, 80), ncol = 3)
expect_equal(toVoxel(fromVoxel(points, image), image), points)

obliqueImage <- denseImage(array(0, c(10L, 10L, 10L)), voxelSize = c(1.2, 0.8, 3), worldTransform = oblique)
expect_equal(toVoxel(fromVoxel(points, obliqueImage), obliqueImage), points)

## Scaled coordinates apply the voxel dimensions but ignore rotation. The
## one-based voxel (3,3,3) is two steps on from the origin voxel (1,1,1), so
## it sits at 2 * voxelSize
expect_equal(as.vector(fromVoxel(c(3, 3, 3), image, type = "scaled")), c(4, 4, 4))
expect_equal(as.vector(toVoxel(c(4, 4, 4), image, type = "scaled")), c(3, 3, 3))

## Voxel coordinates pass through untouched -- already one-based on both sides
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

## --- Geometry objects ------------------------------------------------------

## A bare geometry describes a grid with no data, and every accessor accepts
## one in place of an image
grid <- imageGeometry(c(91L, 109L, 91L), worldTransform = las, unit = "mm")
expect_true(isImageGeometry(grid))
expect_false(isImage(grid))
expect_equal(spatial(grid), 3L)
expect_equal(voxelSize(grid), c(2, 2, 2))
expect_equal(worldTransform(grid), las)
expect_equal(grid@unit, "mm")
expect_identical(geometry(grid), grid)
expect_equal(as.vector(fromVoxel(c(1, 1, 1), grid)), las[1:3, 4])

## Explicit voxel size wins over the one implied by the transform, as for images
expect_equal(voxelSize(imageGeometry(c(4L, 4L, 4L), voxelSize = c(1, 1, 1), worldTransform = las)), c(1, 1, 1))

## Defaults: unit voxels at the origin, unit unknown
plainGrid <- imageGeometry(c(4L, 5L))
expect_equal(voxelSize(plainGrid), c(1, 1))
expect_equal(worldTransform(plainGrid), diag(4))
expect_equal(plainGrid@unit, "unknown")
expect_equal(spatial(imageGeometry()), 0L)

## Validation is the geometry's own
expect_error(imageGeometry(c(4L, 4L), voxelSize = c(1, 1, 1)), "one element per spatial dimension")
expect_error(imageGeometry(c(4L, NA)), "missing or negative")
expect_error(imageGeometry(c(4L, 4L, 4L), worldTransform = sheared), "shear")
expect_error(imageGeometry(4L, unit = c("mm", "cm")), "single value")

## Setters work on a bare geometry, and give back a geometry
voxelSize(grid) <- c(1, 1, 1)
expect_true(isImageGeometry(grid))
expect_equal(worldTransform(grid)[1:3, 4], las[1:3, 4])
worldTransform(grid) <- las
expect_equal(voxelSize(grid), c(2, 2, 2))

## An image's geometry is the same object, and can be replaced wholesale
image <- denseImage(array(0, c(91L, 109L, 91L)), worldTransform = las, unit = "mm")
expect_identical(geometry(image), grid)
expect_identical(geometry(asSparse(image)), grid)
expect_identical(geometry(asPacked(image, "int16")), grid)
geometry(image) <- imageGeometry(c(91L, 109L, 91L))
expect_equal(worldTransform(image), diag(4))
expect_true(isDenseImage(image))
expect_error(geometry(image) <- imageGeometry(c(91L, 109L, 90L)), "does not match")

## A plain array has a default geometry over its leading three dimensions
expect_equal(geometry(array(0, c(3L, 4L, 5L, 6L)))@dims, c(3L, 4L, 5L))
expect_equal(spatial(1:10), 1L)
expect_error(geometry(list(1)), "Cannot find a geometry")

## The world centre of the grid, and its physical extent
expect_equal(centre(grid), as.vector(las %*% c(45, 54, 45, 1))[1:3])
expect_equal(extent(grid), c(182, 218, 182))
expect_equal(centre(imageGeometry(c(5L, 3L))), c(2, 1, 0))
expect_identical(center(grid), centre(grid))
expect_equal(centre(grid), centre(denseImage(array(0, c(91L, 109L, 91L)), geometry = grid)))

## Comparison is by grid and placement, allowing for floating-point noise
expect_true(sameGeometry(grid, image <- denseImage(array(0, c(91L, 109L, 91L)), geometry = grid)))
nudged <- las
nudged[1, 4] <- nudged[1, 4] + 1e-10
expect_true(sameGeometry(grid, imageGeometry(c(91L, 109L, 91L), worldTransform = nudged)))
expect_false(sameGeometry(grid, imageGeometry(c(91L, 109L, 91L), worldTransform = oblique)))
expect_false(sameGeometry(grid, imageGeometry(c(91L, 109L, 90L), worldTransform = las)))

## ...and units only when both are known
expect_true(sameGeometry(grid, imageGeometry(c(91L, 109L, 91L), worldTransform = las)))
expect_false(sameGeometry(grid, imageGeometry(c(91L, 109L, 91L), worldTransform = las, unit = "um")))

## Two-dimensional points are accepted for a two-dimensional grid. Voxel
## coordinates come back with two columns, but world coordinates always have
## three, since a two-dimensional grid may be placed obliquely in world space
flatGrid <- imageGeometry(c(10L, 10L), worldTransform = diag(c(2, 3, 1, 1)))
expect_equal(fromVoxel(c(2, 3), flatGrid), matrix(c(2, 6, 0), 1L))
expect_equal(fromVoxel(c(2, 3), flatGrid, type = "voxel"), matrix(c(2, 3), 1L))
expect_equal(toVoxel(matrix(c(2, 6, 4, 9), 2L, byrow = TRUE), flatGrid), matrix(c(2, 3, 3, 4), 2L, byrow = TRUE))

tilted <- diag(4)
tilted[1:3, 1:3] <- rbind(c(1, 0, 0), c(0, cos(0.5), -sin(0.5)), c(0, sin(0.5), cos(0.5)))
tiltedGrid <- imageGeometry(c(10L, 10L), worldTransform = tilted)
expect_equal(fromVoxel(c(1, 2), tiltedGrid), matrix(c(0, cos(0.5), sin(0.5)), 1L))
expect_equal(toVoxel(fromVoxel(c(4, 7), tiltedGrid), tiltedGrid), matrix(c(4, 7, 1), 1L))

## ...and likewise one-dimensional points for a one-dimensional grid
lineGrid <- imageGeometry(10L, worldTransform = diag(c(2, 1, 1, 1)))
expect_equal(fromVoxel(4, lineGrid), matrix(c(6, 0, 0), 1L))
expect_equal(toVoxel(c(6, 0, 0), lineGrid), matrix(c(4, 1, 1), 1L))
expect_equal(fromVoxel(matrix(c(1, 4), ncol = 1L), lineGrid, type = "voxel"), matrix(c(1, 4), ncol = 1L))

output <- capture.output(print(grid))
expect_true(any(grepl("Image geometry: 91 x 109 x 91", output)))
expect_true(any(grepl("2 x 2 x 2 mm", output)))
