#' World axes and reorientation
#'
#' Every image has an axis order and direction relative to world space.
#' `worldAxes()` reports it as a signed integer for each spatial axis: which
#' world axis it runs along, one to three, negated if its index increases in
#' the negative direction. So `c(-1, 2, 3)` describes an image whose first
#' index increases along the negative first world axis, and whose second and
#' third increase along the second and third. The form is the same as a
#' storage layout's, with world axes in place of storage axes.
#'
#' Nothing here attaches a meaning to the world axes. Conventions that do,
#' such as the anatomical codes used in neuroimaging, belong in packages
#' built on this one.
#'
#' `reorient()` permutes and reverses an image's spatial axes so that they run
#' along the world axes given, updating the geometry to match, so that every
#' voxel keeps its position in world space. Only the 48 axis-aligned
#' orientations are reachable this way; whatever oblique rotation is left over
#' stays in the geometry. For an image with fewer than three spatial
#' dimensions, the entries for the world axes it does not span are ignored.
#'
#' For packed and sparse images this is a change of *view* rather than of
#' data: the image records how its axes map onto its storage, and indexing,
#' `dim()`, the geometry and every apply function work in the view, so no
#' values are moved or copied. A dense image is a plain R array, which base R
#' indexes in memory order, so a dense image is always stored in view order and
#' reorienting one copies its data.
#'
#' @param x An image, or for `worldAxes()` also a geometry.
#' @param to The world axes wanted, in the same form as `worldAxes()` returns:
#'   a signed integer per axis, naming each world axis at most once. The
#'   default puts every axis along its positive world axis, in order.
#' @return `worldAxes()` returns an integer vector with one element per
#'   spatial dimension, up to three. `reorient()` returns an image of the same
#'   class as `x`, with its axes permuted and reversed as needed.
#' @examples
#' image <- denseImage(array(1:24, c(2, 3, 4)), worldTransform = diag(c(-1, 1, 1, 1)))
#' worldAxes(image)
#' worldAxes(reorient(image))
#' @export
worldAxes <- function (x)
{
    geometry <- geometry(x)
    nSpatial <- min(3L, length(geometry@dims))
    if (nSpatial == 0L)
        return(integer(0))

    ## The assignment is the one that aligns the axes best overall, rather
    ## than axis by axis, so that an oblique image near 45 degrees cannot have
    ## two voxel axes claim the same world axis
    block <- geometry@orientation[1:3, seq_len(nSpatial), drop = FALSE]
    candidates <- as.matrix(expand.grid(rep(list(1:3), nSpatial)))
    candidates <- candidates[apply(candidates, 1L, anyDuplicated) == 0L, , drop = FALSE]
    scores <- apply(candidates, 1L, function (world) sum(abs(block[cbind(world, seq_len(nSpatial))])))
    world <- as.integer(candidates[which.max(scores), ])
    ifelse(block[cbind(world, seq_len(nSpatial))] >= 0, world, -world)
}

#' @rdname worldAxes
#' @export
reorient <- function (x, to = 1:3)
{
    if (!isImage(x))
        stop("Only an image can be reoriented")

    to <- checkWorldAxes(to)
    geometry <- geometry(x)
    nSpatial <- min(3L, length(geometry@dims))
    axes <- worldAxes(geometry)

    ## The image's axes are put in the order their world axes take in the
    ## target, and any pointing the other way are reversed
    if (!all(abs(axes) %in% abs(to)))
        stop("The target does not say which way every spatial axis should run")
    order <- order(match(abs(axes), abs(to)))
    reversed <- sign(axes[order]) != sign(to[match(abs(axes[order]), abs(to))])

    if (identical(order, seq_len(nSpatial)) && !any(reversed))
        return(x)

    nDims <- length(dim(x))
    permutation <- c(order, seq_len(nDims)[-seq_len(nSpatial)])
    flips <- c(reversed, rep(FALSE, nDims - nSpatial))
    newGeometry <- permuteGeometry(geometry, order, reversed)

    if (isDenseImage(x))
    {
        ## Relative to its current order, which is the storage order of a
        ## dense image, the new view is simply the permutation with signs
        layout <- ifelse(flips, -permutation, permutation)
        values <- viewGather(as.array(x), dim(x)[permutation], layout)
        return(denseImage(values, geometry = newGeometry))
    }

    ## A view axis maps to whatever storage axis the old one it came from did,
    ## and reversing it reverses that
    layout <- x@layout[permutation] * ifelse(flips, -1L, 1L)
    S7::set_props(x, dims = x@dims[permutation], layout = as.integer(layout), geometry = newGeometry)
}

checkWorldAxes <- function (to)
{
    if (!is.numeric(to) || length(to) < 1L || length(to) > 3L || anyNA(to) || any(to != round(to)))
        stop("World axes must be given as one to three signed integers, such as c(1, 2, 3)")
    to <- as.integer(to)
    if (any(to == 0L) || any(abs(to) > 3L))
        stop("World axes must be numbered from 1 to 3, negated for the negative direction")
    if (anyDuplicated(abs(to)))
        stop("The target names the same world axis twice")
    to
}

## Permutes and reverses the spatial axes of a geometry, keeping every voxel
## where it is in world space. Reversing an axis moves the origin to what was
## its last voxel. The orientation's columns for any axes beyond the spatial
## ones follow the retained ones, keeping the block a rotation or reflection
permuteGeometry <- function (geometry, order, reversed)
{
    n <- length(order)
    xform <- worldTransform(geometry)
    orientation <- geometry@orientation

    for (i in which(reversed))
    {
        axis <- order[i]
        orientation[1:3, 4] <- orientation[1:3, 4] + xform[1:3, axis] * (geometry@dims[axis] - 1)
    }

    columns <- c(order, setdiff(1:3, order))
    block <- orientation[1:3, columns, drop = FALSE]
    block[, seq_len(n)] <- block[, seq_len(n), drop = FALSE] %*% diag(ifelse(reversed, -1, 1), n)
    orientation[1:3, 1:3] <- block

    rest <- seq_along(geometry@dims)[-seq_len(n)]
    S7::set_props(geometry, dims = geometry@dims[c(order, rest)],
                  voxelSize = geometry@voxelSize[c(order, rest)], orientation = orientation)
}

## Layouts. A layout gives, for each view axis, the one-based rank of the
## storage axis it corresponds to, negated where the view runs backwards along
## it; see View.h for the full description

isIdentityLayout <- function (layout) identical(as.integer(layout), seq_along(layout))

## What to pass to compiled code: nothing at all for the identity, so that the
## common case takes the direct path
layoutArg <- function (x) if (isIdentityLayout(x@layout)) NULL else x@layout

## One-based indices in the image's own order, as indices into its storage
storageIndices <- function (x, i)
{
    i <- as.double(i)
    if (isIdentityLayout(x@layout)) i else viewIndices(i, x@dims, x@layout)
}

printLayout <- function (layout)
{
    if (!isIdentityLayout(layout))
        cat(sprintf("  Storage layout     : %s\n", paste(layout, collapse = ", ")))
}

checkLayout <- function (layout, n)
{
    if (length(layout) != n)
        return(paste0("@layout must have one entry per dimension (", n, ")"))
    if (anyNA(layout) || any(layout == 0L) || !setequal(abs(layout), seq_len(n)))
        return(paste0("@layout must name each storage axis from 1 to ", n, " exactly once"))
    NULL
}

#' @rdname worldAxes
#' @details `storageLayout()` returns the layout of an image: one entry per
#'   axis, giving the one-based rank of the storage axis it corresponds to
#'   (rank one varying fastest), negated where the image runs along that axis
#'   backwards. It is the identity for a dense image, and for any image that
#'   has not been reoriented or read from a file stored in another order.
#' @export
storageLayout <- function (x)
{
    if (isPackedImage(x) || isSparseImage(x))
        x@layout
    else
        seq_along(dim(x) %||% length(x))
}
