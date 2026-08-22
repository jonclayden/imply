#' Image geometry
#'
#' Accessors for the geometry of an image: the number of spatial dimensions,
#' the voxel size and the placement of the image in world space.
#'
#' These are deliberately free of any dependency on a file format, and of any
#' assumption that geometry is anatomical: an image with no meaningful spatial
#' interpretation is free to leave it at the identity default and never think
#' about it again.
#'
#' Voxel size and world placement are stored, and can be set, independently of
#' one another, unlike a NIfTI xform, which bakes voxel size into the same
#' matrix that carries rotation and translation and so has to be kept in sync
#' by hand. Here, `voxelSize<-` never touches rotation or translation, and
#' the stored placement is always a rigid transform (rotation or reflection
#' plus translation; no scale and no shear) so that the two cannot drift out
#' of agreement.
#'
#' `worldTransform()` composes the two into the single 4x4 affine that other
#' packages expect. Setting it back decomposes the matrix into rotation and
#' scale; a matrix that doesn't decompose that way (i.e. one with genuine
#' shear) is rejected rather than silently mangled. In practice a sheared
#' NIfTI `sform` usually means the field is being used to carry an affine
#' registration or normalisation result (to Talairach or MNI space, say)
#' rather than to describe voxel storage geometry, which is a different kind
#' of information than this package models.
#'
#' `toVoxel()` and `fromVoxel()` use one-based voxel coordinates, matching
#' `x[i, j, k]` indexing. The affine itself, and the rest of the package's C++
#' API, is zero-based throughout.
#'
#' @param x An image, or for `worldTransform` either an image or a 4x4 matrix.
#' @param value A replacement value.
#' @param points A matrix of points, one per row and three columns. Where
#'   `type` is `"voxel"` (the default convention for `fromVoxel()`'s input and
#'   `toVoxel()`'s output), these are one-based, as for `x[i, j, k]`.
#' @param type The coordinate convention of `points`: `"voxel"` (one-based),
#'   `"scaled"` (millimetres from the one-based origin, ignoring rotation) or
#'   `"world"` (fully transformed).
#' @param round Rounding strategy: `"none"`, `"conventional"` for nearest
#'   neighbour, or `"probabilistic"` for a stochastic nearest neighbour with
#'   probability proportional to proximity.
#' @param bounds Optional image extents, used only by probabilistic rounding to
#'   avoid selecting a location beyond the end of the image.
#' @name geometry
NULL

`%||%` <- function (x, y) if (is.null(x)) y else x

#' @rdname geometry
#' @export
isImage <- function (x) isDenseImage(x) || isSparseImage(x) || isPackedImage(x)

#' @rdname geometry
#' @export
spatial <- function (x) attr(x, "spatial") %||% min(3L, length(dim(x)))

#' @rdname geometry
#' @export
voxelSize <- function (x) attr(x, "voxelSize") %||% rep(1, spatial(x))

#' @rdname geometry
#' @export
`voxelSize<-` <- function (x, value)
{
    value <- as.double(value)
    if (!isImage(x))
        x <- denseImage(x)
    x@voxelSize <- value
    x
}

#' @rdname geometry
#' @export
worldTransform <- function (x)
{
    ## isImage() is checked first, and matrix-ness second, because a
    ## two-dimensional image is itself a matrix: without this order such an
    ## image would be misread as a raw transform to validate rather than an
    ## image whose transform is wanted
    if (isImage(x))
        composeTransform(attr(x, "orientation") %||% diag(4), voxelSize(x))
    else if (is.matrix(x))
        validateXform(x)
    else
        composeTransform(diag(4), voxelSize(x))
}

#' @rdname geometry
#' @export
`worldTransform<-` <- function (x, value)
{
    ## As for voxelSize<-(): preserve whatever image class x already is
    if (!isImage(x))
        x <- denseImage(x)
    decomposed <- decomposeTransform(validateXform(value), x@spatial)
    x@orientation <- decomposed$orientation
    x@voxelSize <- decomposed$voxelSize
    x
}

validateXform <- function (value)
{
    value <- as.matrix(value)
    storage.mode(value) <- "double"
    if (!identical(dim(value), c(4L, 4L)))
        stop("Transform must be a 4x4 matrix")
    if (anyNA(value))
        stop("Transform must not contain missing values")
    if (!isTRUE(all.equal(value[4, ], c(0, 0, 0, 1))))
        stop("Transform must be affine, with a final row of (0, 0, 0, 1)")
    dimnames(value) <- NULL
    value
}

## The tolerance below is deliberately loose relative to floating-point noise:
## NIfTI sform/qform fields are stored as float32 in the header, which can
## introduce relative error of order 1e-7 in each entry, and considerably more
## after normalising a column by a voxel size of only a few millimetres. A
## genuine shear, in contrast, comes from something like a 12-parameter affine
## registration and is essentially never this close to orthogonal
orthogonalityTolerance <- 1e-4

## Splits a full affine into a rigid placement (unit-length, mutually
## orthogonal columns, i.e. a rotation or reflection, plus translation) and a
## voxel size. Errors if the 3x3 block doesn't decompose that way, which in
## practice usually means it encodes a general affine registration rather than
## voxel storage geometry - a different kind of information than an
## orientation matrix here is meant to hold
decomposeTransform <- function (xform, spatial)
{
    block <- xform[1:3, 1:3, drop = FALSE]
    norms <- sqrt(colSums(block^2))
    if (any(norms[seq_len(spatial)] < .Machine$double.eps^0.5))
        stop("Transform has a zero-length axis and cannot be decomposed")

    unit <- block
    for (i in 1:3)
        unit[, i] <- if (norms[i] > 0) block[, i] / norms[i] else block[, i]

    crossTerms <- crossprod(unit) - diag(3)
    if (max(abs(crossTerms[upper.tri(crossTerms)])) > orthogonalityTolerance)
        stop("Transform cannot be decomposed into rotation and voxel size ",
             "(it contains shear); resample the image, or supply a rigid ",
             "transform with voxelSize<-/worldTransform<- set separately")

    ## Axes beyond `spatial` have nowhere to store a non-unit scale
    if (spatial < 3L && any(abs(norms[(spatial + 1L):3] - 1) > orthogonalityTolerance))
        stop("Transform implies a non-unit scale on an axis beyond the ",
             "image's spatial dimensions, which cannot be represented")

    orientation <- diag(4)
    orientation[1:3, 1:3] <- unit
    orientation[1:3, 4] <- xform[1:3, 4]

    list(orientation = orientation, voxelSize = norms[seq_len(spatial)])
}

## The inverse of decomposeTransform(): recombines a rigid placement and a
## voxel size into the single affine other packages expect
composeTransform <- function (orientation, voxelSize)
{
    result <- orientation
    n <- min(3L, length(voxelSize))
    if (n > 0L)
    {
        for (i in seq_len(n))
            result[1:3, i] <- result[1:3, i] * voxelSize[i]
    }
    result
}

#' @rdname geometry
#' @export
toVoxel <- function (points, x, type = "world", round = "none", bounds = NULL)
{
    points <- asPointMatrix(points)
    result <- pointsToVoxel(points, worldTransform(x), voxelSize(x), type)
    if (!identical(round, "none"))
    {
        if (is.null(bounds) && !is.matrix(x))
            bounds <- as.double(dim(x)[seq_len(min(3L, spatial(x)))])
        result <- roundPoints(result, round, bounds)
    }
    ## Standard one-based indexing for R
    if (!identical(type, "voxel"))
        result <- result + 1
    result
}

#' @rdname geometry
#' @export
fromVoxel <- function (points, x, type = "world")
{
    points <- asPointMatrix(points)
    ## The inverse of toVoxel()'s shift: points arrive one-based and are
    ## brought back to zero-based before reaching the affine, unless they are
    ## staying in voxel space, in which case there is nothing to shift
    if (!identical(type, "voxel"))
        points <- points - 1
    pointsFromVoxel(points, worldTransform(x), voxelSize(x), type)
}

asPointMatrix <- function (points)
{
    if (is.null(dim(points)))
        points <- matrix(points, nrow = 1L)
    points <- as.matrix(points)
    storage.mode(points) <- "double"
    if (ncol(points) != 3L)
        stop("Points must be given as a matrix with three columns")
    points
}
