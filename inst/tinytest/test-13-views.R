## Orientation, reorientation and views. The reference throughout is the dense
## image: a packed or sparse image whose axes have been permuted and reversed
## by changing its view must behave exactly like a dense image whose data were
## actually moved.

set.seed(13)
las <- rbind(c(-2, 0, 0, 90), c(0, 2.5, 0, -126), c(0, 0, 3, -72), c(0, 0, 0, 1))
values <- array(round(rnorm(4 * 5 * 6 * 3), 2), c(4L, 5L, 6L, 3L))
values[abs(values) < 0.8] <- 0
dense <- denseImage(values, worldTransform = las, unit = "mm")

## --- World axes ------------------------------------------------------------

## One signed world axis per spatial axis: here the first index increases
## along the negative first world axis
expect_identical(worldAxes(dense), c(-1L, 2L, 3L))
expect_identical(worldAxes(geometry(dense)), c(-1L, 2L, 3L))
expect_identical(worldAxes(denseImage(array(0, c(2, 2, 2)))), 1:3)
expect_identical(worldAxes(denseImage(matrix(0, 2, 3), worldTransform = diag(c(1, -1, 1, 1)))), c(1L, -2L))
expect_identical(worldAxes(imageGeometry()), integer(0))

## A permuted grid reports which world axis each voxel axis follows
swapped <- diag(4)
swapped[1:3, 1:3] <- rbind(c(0, 0, 1), c(-1, 0, 0), c(0, 1, 0))
expect_identical(worldAxes(imageGeometry(c(2L, 2L, 2L), worldTransform = swapped)), c(-2L, 3L, 1L))

## Near 45 degrees the axes are assigned jointly, so two voxel axes can never
## claim the same world axis
tilt <- diag(4)
angle <- pi / 4 - 0.01
tilt[1:2, 1:2] <- rbind(c(cos(angle), -sin(angle)), c(sin(angle), cos(angle)))
expect_equal(sort(abs(worldAxes(imageGeometry(c(3L, 3L, 3L), worldTransform = tilt)))), 1:3)

expect_error(reorient(dense, c(1, 4, 2)), "numbered from 1 to 3")
expect_error(reorient(dense, c(1, 0, 2)), "numbered from 1 to 3")
expect_error(reorient(dense, c(1, -1, 3)), "same world axis twice")
expect_error(reorient(dense, c(1.5, 2, 3)), "signed integers")
expect_error(reorient(dense, "RAS"), "signed integers")
expect_error(reorient(dense, c(1, 2)), "every spatial axis")
expect_error(reorient(geometry(dense)), "Only an image")

## --- Reorientation against the dense reference -----------------------------

## World position of a one-based voxel of the original
worldOf <- function (voxel) as.vector(las %*% c(voxel - 1, 1))[1:3]

targets <- list(c(1L, 2L, 3L), c(-1L, -2L, -3L), c(3L, 2L, 1L), c(-3L, -1L, -2L), c(-1L, 2L, 3L))
for (to in targets)
{
    info <- paste("reorient to", paste(to, collapse = ", "))
    reference <- reorient(dense, to)
    packed <- reorient(asPacked(dense, "float32"), to)
    sparse <- reorient(asSparse(dense), to)

    expect_true(isDenseImage(reference), info = info)
    expect_true(isPackedImage(packed), info = info)
    expect_true(isSparseImage(sparse), info = info)

    for (image in list(reference, packed, sparse))
    {
        expect_identical(worldAxes(image), to, info = info)
        expect_equal(centre(image), centre(dense), info = info)
        expect_equal(geometry(image)@unit, "mm", info = info)
    }

    ## Every voxel keeps its place in world space, and its value
    for (voxel in list(c(1, 1, 1), c(2, 3, 4), c(4, 5, 6)))
    {
        moved <- round(toVoxel(worldOf(voxel), reference))
        expect_equal(as.array(reference)[cbind(moved, 2)], values[cbind(t(voxel), 2)], info = info)
    }

    ## Views read exactly what the moved data hold
    expect_identical(as.array(sparse), as.array(reference), info = info)
    expect_equal(as.array(packed), as.array(reference), tolerance = 1e-6, info = info)
    expect_identical(dim(packed), dim(reference), info = info)

    ## Element access maps through the view without materialising
    picks <- c(1, 7, 55, 200, 359)
    expect_identical(sparse[picks], as.array(reference)[picks], info = info)
    expect_equal(packed[picks], as.array(reference)[picks], tolerance = 1e-6, info = info)
    expect_identical(sparse[cbind(2, 3, 1, 3)], as.array(reference)[cbind(2, 3, 1, 3)], info = info)

    ## The apply and reduce engines walk the view
    expect_identical(as.array(voxelApply(sparse, sum)), as.array(voxelApply(reference, sum)), info = info)
    expect_identical(imapply(sparse, c(1, 3), max), imapply(as.array(reference), c(1, 3), max), info = info)
    expect_identical(imapply(sparse, 4, sum), imapply(as.array(reference), 4, sum), info = info)
    expect_equal(imreduce(packed, 2:3, "sum"), imreduce(as.array(packed), 2:3, "sum"), info = info)
    expect_identical(imreduce(sparse, 1, "which.max"), imreduce(as.array(reference), 1, "which.max"), info = info)
    expect_identical(as.array(lineApply(sparse, sum, axis = 2)), as.array(lineApply(reference, sum, axis = 2)), info = info)

    ## The mask and the masked matrix are in view order too
    expect_identical(mask(sparse), apply(as.array(reference) != 0, 1:3, any), info = info)
    expect_identical(maskedMatrix(sparse), maskedMatrix(asSparse(reference)), info = info)
    stored <- storedValues(sparse)
    expect_identical(stored$values[, order(stored$locations), drop = FALSE], maskedMatrix(sparse), info = info)
    expect_identical(sort(stored$locations), which(mask(sparse)), info = info)

    ## Reorienting back restores the original exactly
    expect_identical(as.array(reorient(sparse, c(-1L, 2L, 3L))), values, info = info)
    expect_identical(as.array(reorient(reference, c(-1L, 2L, 3L))), values, info = info)
    expect_equal(worldTransform(reorient(packed, c(-1L, 2L, 3L))), las, info = info)
}

## A view never moves data: the stored bytes are the same object's
packed <- asPacked(dense, "int16")
turned <- reorient(packed, c(-2L, -3L, 1L))
expect_identical(turned@values, packed@values)
expect_false(identical(storageLayout(turned), 1:4))
expect_identical(storageLayout(dense), 1:4)
expect_identical(storageLayout(reorient(dense, c(-2L, -3L, 1L))), 1:4)

## An image already in the orientation asked for comes back unchanged
expect_identical(reorient(dense, c(-1L, 2L, 3L)), dense)

## Values held at each location travel with it, and the value dimensions stay
## put
expect_equal(dim(reorient(dense, c(3L, 2L, 1L))), c(6L, 5L, 4L, 3L))

## Two-dimensional images use the entries for the world axes they span, and
## the default puts every axis along its positive world axis
flat <- denseImage(matrix(1:6, 2, 3), worldTransform = diag(c(1, -1, 1, 1)))
expect_identical(worldAxes(reorient(flat)), 1:2)
expect_identical(as.array(reorient(flat)), matrix(1:6, 2, 3)[, 3:1])
expect_identical(worldAxes(reorient(flat, c(-2L, -3L, -1L))), c(-2L, -1L))

## --- Views in other operations ---------------------------------------------

sparse <- reorient(asSparse(dense), c(3L, -2L, 1L))
reference <- reorient(dense, c(3L, -2L, 1L))

## Arithmetic on the stored values keeps the view
doubled <- sparse * 2
expect_true(isSparseImage(doubled))
expect_identical(storageLayout(doubled), storageLayout(sparse))
expect_identical(as.array(doubled), as.array(reference) * 2)

## Two sparse images stored in different orders are combined densely, in
## view order, rather than location by location
other <- reorient(asSparse(reorient(dense, c(3L, -2L, 1L))), c(3L, -2L, 1L))
expect_false(identical(storageLayout(other), storageLayout(sparse)))
expect_equal(as.array(sparse + other), as.array(reference) * 2)
expect_equal(as.array(sparse * other), as.array(reference)^2)

## Converting keeps the view's order
expect_identical(as.array(asDense(sparse)), as.array(reference))
expect_identical(as.array(asSparse(asDense(sparse))), as.array(reference))
expect_equal(as.array(asPacked(sparse, "float32")), as.array(reference), tolerance = 1e-6)

## Masked application sees the same locations in the same order
brain <- mask(sparse)
expect_equal(as.array(voxelApply(sparse, mean, mask = brain)), as.array(voxelApply(reference, mean, mask = brain)))

## --- Validation ------------------------------------------------------------

expect_error(packedImage(packed@values, "int16", dims = dim(packed), layout = c(1L, 1L, 2L, 3L)), "exactly once")
expect_error(packedImage(packed@values, "int16", dims = dim(packed), layout = 1:3), "one entry per dimension")
expect_error(sparseImage(mask = asSparse(dense)@mask, values = asSparse(dense)@values, dim = dim(dense),
                         layout = c(4L, 2L, 3L, 1L)), "only reorder its spatial axes")

output <- capture.output(print(turned))
expect_true(any(grepl("Storage layout", output)))
expect_false(any(grepl("Storage layout", capture.output(print(packed)))))
