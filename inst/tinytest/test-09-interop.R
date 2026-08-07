## Bridges to the image classes of packages that handle file formats, and the
## public C++ header interface.

las <- rbind(c(-2, 0, 0,   90),
             c( 0, 2, 0, -126),
             c( 0, 0, 2,  -72),
             c( 0, 0, 0,    1))

set.seed(1)
values <- array(rnorm(6L * 7L * 8L), c(6L, 7L, 8L))
image <- denseImage(values, pixdim = c(2, 2, 2), xform = las)

## --- RNifti ----------------------------------------------------------------

## A niftiImage is an array carrying RNifti's own class and attributes, so
## comparisons are made against the bare values
bareValues <- function (x) {
    x <- as.array(x)
    attributes(x) <- list(dim = dim(x))
    x
}

if (requireNamespace("RNifti", quietly = TRUE))
{
    nifti <- toNifti(image)
    expect_true(inherits(nifti, "niftiImage"))

    ## A NIfTI xform maps voxel to world coordinates, which is exactly what
    ## imply holds, so nothing has to be reinterpreted
    expect_equal(unname(RNifti::xform(nifti)[1:4, 1:4]), las)
    expect_equal(RNifti::pixdim(nifti)[1:3], c(2, 2, 2))
    expect_equal(RNifti::orientation(nifti), orientation(image))
    expect_equal(bareValues(nifti), values)

    back <- fromNifti(nifti)
    expect_true(isDenseImage(back))
    expect_equal(as.array(back), values)
    expect_equal(xform(back), las)
    expect_equal(pixdim(back), c(2, 2, 2))
    expect_equal(orientation(back), "LAS")
    expect_equal(spatial(back), 3L)

    ## Anisotropic voxels, where a transform and a pixdim that disagree would
    ## show up
    oblique <- rbind(c(-1, 0, 0, 7), c(0, 3, 0, -2), c(0, 0, 5, 4), c(0, 0, 0, 1))
    anisotropic <- denseImage(array(rnorm(64), c(4L, 4L, 4L)), pixdim = c(1, 3, 5), xform = oblique)
    converted <- toNifti(anisotropic)
    expect_equal(unname(RNifti::xform(converted)[1:4, 1:4]), oblique)
    expect_equal(RNifti::pixdim(converted)[1:3], c(1, 3, 5))
    expect_equal(xform(fromNifti(converted)), oblique)
    expect_equal(pixdim(fromNifti(converted)), c(1, 3, 5))

    ## An image RNifti made itself, rather than one round-tripped from here
    native <- RNifti::asNifti(array(rnorm(64), c(4L, 4L, 4L)))
    expect_true(isDenseImage(fromNifti(native)))
    expect_equal(dim(fromNifti(native)), c(4L, 4L, 4L))
    expect_equal(orientation(fromNifti(native)), RNifti::orientation(native))

    ## Nothing of RNifti's bookkeeping survives into the converted image
    expect_null(attr(as.array(fromNifti(nifti)), "pixunits"))
    expect_null(attr(as.array(fromNifti(nifti)), ".nifti_image_ptr"))
    expect_false(inherits(fromNifti(nifti), "niftiImage"))

    ## Sparse and packed images convert by materialising first
    expect_equal(bareValues(toNifti(asSparse(image))), values)
    expect_equal(bareValues(toNifti(asPacked(image, "float32"))), values, tolerance = 1e-6)
}

## --- tractor.base ----------------------------------------------------------

if (requireNamespace("tractor.base", quietly = TRUE))
{
    mri <- toMriImage(image)
    expect_true(inherits(mri, "MriImage"))
    expect_equal(as.array(mri), values)
    expect_equal(abs(mri$getVoxelDimensions()), c(2, 2, 2))

    recovered <- fromMriImage(mri)
    expect_true(isDenseImage(recovered))
    expect_equal(as.array(recovered), values)
    expect_equal(pixdim(recovered), c(2, 2, 2))
    expect_equal(dim(recovered), c(6L, 7L, 8L))
}

## --- The public C++ interface ----------------------------------------------

## Everything is header-only, so a package linking to imply needs no library.
## These checks are on the shipped layout rather than on behaviour
headers <- system.file("include", package = "imply")
expect_true(nzchar(headers))
expect_true(file.exists(file.path(headers, "imply.h")))

for (header in c("Raster.h", "Space.h", "Blocks.h", "Parallel.h", "Storage.h",
                 "Dispatch.h", "RImage.h", "Sparse.h", "Narrow.h", "Sink.h"))
    expect_true(file.exists(file.path(headers, "imply", header)), info = header)

## The umbrella pulls in every one of them
umbrella <- readLines(file.path(headers, "imply.h"))
for (header in c("Raster.h", "Space.h", "Blocks.h", "Parallel.h", "Sparse.h", "Narrow.h"))
    expect_true(any(grepl(paste0('include "imply/', header), umbrella, fixed = TRUE)), info = header)

## Space.h carries its own definitions, so there is nothing to link against
space <- readLines(file.path(headers, "imply", "Space.h"))
expect_true(any(grepl("^inline .*ImageSpace::orientation", space)))
expect_true(any(grepl("^inline .*Affine::inverse", space)))

## Raster.h is free of any dependency on R, so it can be used on its own
raster <- readLines(file.path(headers, "imply", "Raster.h"))
expect_false(any(grepl("Rcpp\\.h|Rinternals\\.h|<R\\.h>", raster)))
