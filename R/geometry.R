#' Image geometry
#'
#' An image geometry describes the spatial grid an image is sampled on,
#' without any of its data: the extent of the grid, the voxel size, the
#' placement of the grid in world space and the unit of measurement. Every
#' image holds one, and any of the functions here may be given either an image
#' or a bare geometry. A bare geometry is useful in its own right, to describe
#' a space that no image has been created in yet, such as the target grid of a
#' resampling.
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
#' @param dims The extent of the grid, one per spatial dimension.
#' @param voxelSize Voxel size, one per spatial dimension.
#' @param worldTransform A 4x4 affine transform mapping (zero-based) voxel to
#'   world coordinates, from which the voxel size is also taken unless it is
#'   given separately.
#' @param unit The unit of measurement of voxel size and world coordinates.
#' @param x,y Images or geometries. Anything else with a dimension attribute is
#'   treated as having unit voxels at the origin, with the leading three
#'   dimensions (or all of them, if fewer) spatial. `worldTransform()` also
#'   accepts a bare 4x4 matrix, which is validated and returned.
#' @param value A replacement value.
#' @param tolerance Numerical tolerance for comparing world transforms.
#' @param points A matrix of points, one per row and one to three columns.
#'   Missing columns are taken to be zero. Where `type` is `"voxel"` (the
#'   default convention for `fromVoxel()`'s input and `toVoxel()`'s output),
#'   these are one-based, as for `x[i, j, k]`.
#' @param type The coordinate convention of `points`: `"voxel"` (one-based),
#'   `"scaled"` (millimetres from the one-based origin, ignoring rotation) or
#'   `"world"` (fully transformed).
#' @param round Rounding strategy: `"none"`, `"conventional"` for nearest
#'   neighbour, or `"probabilistic"` for a stochastic nearest neighbour with
#'   probability proportional to proximity.
#' @param bounds Optional image extents, used only by probabilistic rounding to
#'   avoid selecting a location beyond the end of the image.
#' @return `imageGeometry()` and `geometry()` return an object of S7 class
#'   `imageGeometry`, with properties `dims`, `voxelSize`, `orientation` (the
#'   rigid part of the world transform) and `unit`. `isImage()` and
#'   `isImageGeometry()` return a Boolean value indicating whether their
#'   argument is one of the package's image types, or a geometry,
#'   respectively. `spatial()` returns the number of spatial dimensions.
#'   `voxelSize()` returns a vector of sizes in each spatial dimension.
#'   `worldTransform()` returns a numeric affine transform matrix. `centre()`
#'   returns the world coordinates of the centre of the grid, always as three
#'   values; `center()` is an alias. `extent()` returns the physical size of
#'   the grid along each spatial axis. `sameGeometry()` returns a Boolean
#'   value, which is `TRUE` if the two grids are the same size and are placed
#'   identically in world space, up to `tolerance`; units are compared only if
#'   both are known.
#'   `toVoxel()` and `fromVoxel()` return matrices of transformed points, one
#'   per row. Voxel coordinates have as many columns as `points`, and world
#'   or scaled coordinates always have three. The assignment functions are
#'   called for their side-effects.
#' @name geometry
NULL

`%||%` <- function (x, y) if (is.null(x)) y else x

#' @rdname geometry
#' @export
imageGeometry <- S7::new_class("imageGeometry",
    properties = list(
        dims = S7::class_integer,
        voxelSize = S7::class_double,
        orientation = S7::class_double,
        unit = S7::class_character
    ),
    validator = function (self) {
        if (anyNA(self@dims) || any(self@dims < 0L))
            return("@dims must not be missing or negative")

        if (length(self@voxelSize) != length(self@dims))
            return(paste0("@voxelSize must have one element per spatial dimension (", length(self@dims), ")"))
        if (anyNA(self@voxelSize))
            return("@voxelSize must not be missing")
        if (any(self@voxelSize <= 0))
            return("@voxelSize must be strictly positive")

        if (!identical(dim(self@orientation), c(4L, 4L)))
            return("@orientation must be a 4x4 matrix")
        if (anyNA(self@orientation))
            return("@orientation must not contain missing values")
        if (!isTRUE(all.equal(self@orientation[4, ], c(0, 0, 0, 1))))
            return("@orientation must be affine, with a final row of (0, 0, 0, 1)")
        block <- self@orientation[1:3, 1:3, drop = FALSE]
        if (max(abs(crossprod(block) - diag(3))) > orthogonalityTolerance)
            return("@orientation must be rigid: a rotation or reflection, with no scale or shear")

        if (length(self@unit) != 1L || is.na(self@unit))
            return("@unit must be a single value")

        NULL
    },
    constructor = function (dims = integer(0), voxelSize = NULL, worldTransform = NULL, unit = NULL)
    {
        dims <- as.integer(dims)
        decomposed <- if (is.null(worldTransform)) NULL
                      else decomposeTransform(validateXform(worldTransform), length(dims))

        S7::new_object(S7::S7_object(),
            dims = dims,
            voxelSize = as.double(voxelSize %||% decomposed$voxelSize %||% rep(1, length(dims))),
            orientation = decomposed$orientation %||% diag(4),
            unit = as.character(unit %||% "unknown"))
    })

S7::S4_register(imageGeometry)

S7::method(print, imageGeometry) <- function (x, ...)
{
    cat(sprintf("Image geometry: %s\n", if (length(x@dims) == 0L) "no spatial dimensions"
                                        else paste(x@dims, collapse = " x ")))
    printGeometry(x)
    invisible(x)
}

## The lines describing a grid, shared by the print methods of every image
printGeometry <- function (geometry)
{
    if (length(geometry@dims) > 0L)
    {
        cat(sprintf("  Spatial dimensions : %s\n", paste(geometry@dims, collapse = " x ")))
        cat(sprintf("  Voxel size         : %s %s\n", paste(signif(geometry@voxelSize, 4), collapse = " x "),
                    if (geometry@unit == "unknown") "(unit unknown)" else geometry@unit))
    }
}

#' @rdname geometry
#' @export
isImageGeometry <- function (x) S7::S7_inherits(x, imageGeometry)

#' @rdname geometry
#' @export
isImage <- function (x) isDenseImage(x) || isSparseImage(x) || isPackedImage(x)

#' @rdname geometry
#' @export
geometry <- function (x)
{
    if (isImageGeometry(x))
        return(x)
    if (isImage(x))
        return(x@geometry)
    if (!is.atomic(x))
        stop("Cannot find a geometry for an object of class ", class(x)[1L])

    dims <- dim(x) %||% length(x)
    imageGeometry(dims[seq_len(min(3L, length(dims)))])
}

#' @rdname geometry
#' @export
`geometry<-` <- function (x, value)
{
    value <- geometry(value)
    if (isImageGeometry(x))
        return(value)
    if (!isImage(x))
        x <- denseImage(x)
    x@geometry <- value
    x
}

## Replaces some properties of a geometry, or of the geometry an image holds,
## all at once, so that validation sees only the final state. Whatever image
## class x already is is preserved, so setters never densify as a side effect;
## anything that is not yet an image becomes a dense one
updateGeometry <- function (x, ...)
{
    if (isImageGeometry(x))
        return(S7::set_props(x, ...))
    if (!isImage(x))
        x <- denseImage(x)
    x@geometry <- S7::set_props(x@geometry, ...)
    x
}

## Builds the geometry of an image under construction. Explicit arguments win,
## then a value implied by another explicit argument (a worldTransform implies
## both orientation and voxel size), then the geometry of `from`, then defaults
resolveGeometry <- function (dims, spatial = NULL, voxelSize = NULL, worldTransform = NULL,
                             unit = NULL, from = NULL)
{
    nDims <- length(dims)
    if (!is.null(from))
        from <- geometry(from)

    spatial <- as.integer(spatial %||% if (is.null(from)) min(3L, nDims) else length(from@dims))
    if (length(spatial) != 1L || is.na(spatial))
        stop("The number of spatial dimensions must be a single value")
    if (spatial < 0L || spatial > nDims)
        stop("The number of spatial dimensions must be between 0 and ", nDims)
    spatialDims <- as.integer(dims[seq_len(spatial)])

    if (is.null(from))
        return(imageGeometry(spatialDims, voxelSize, worldTransform, unit))

    if (!identical(from@dims, spatialDims))
        stop("The geometry given is for a grid of ", formatDims(from@dims),
             ", but the image's spatial dimensions are ", formatDims(spatialDims))

    changes <- list()
    if (!is.null(worldTransform))
        changes <- decomposeTransform(validateXform(worldTransform), spatial)
    if (!is.null(voxelSize))
        changes$voxelSize <- as.double(voxelSize)
    if (!is.null(unit))
        changes$unit <- as.character(unit)
    do.call(S7::set_props, c(list(from), changes))
}

## Removes spatial axes from a geometry, keeping the rest in order. A location
## in the result stands for a whole line or plane across the dropped axes, so
## it is placed at the centre of that line or plane. The columns of the rigid
## block are reordered rather than removed, with the dropped axes filling the
## slots beyond the retained ones, so that the block stays a rotation or
## reflection and those axes keep the unit scale decomposeTransform() expects
dropAxes <- function (geometry, axes)
{
    n <- length(geometry@dims)
    keep <- setdiff(seq_len(n), axes)

    xform <- worldTransform(geometry)
    orientation <- geometry@orientation
    for (axis in axes[axes <= 3L])
        orientation[1:3, 4] <- orientation[1:3, 4] + xform[1:3, axis] * (geometry@dims[axis] - 1) / 2

    columns <- c(keep, axes, setdiff(1:3, c(keep, axes)))
    orientation[1:3, 1:3] <- orientation[1:3, columns[columns <= 3L], drop = FALSE]

    S7::set_props(geometry, dims = geometry@dims[keep], voxelSize = geometry@voxelSize[keep],
                  orientation = orientation)
}

formatDims <- function (dims)
    if (length(dims) == 0L) "no dimensions" else paste(dims, collapse = " x ")

#' @rdname geometry
#' @export
spatial <- function (x) length(geometry(x)@dims)

#' @rdname geometry
#' @export
voxelSize <- function (x) geometry(x)@voxelSize

#' @rdname geometry
#' @export
`voxelSize<-` <- function (x, value) updateGeometry(x, voxelSize = as.double(value))

#' @rdname geometry
#' @export
worldTransform <- function (x)
{
    ## Images and geometries are checked first, and matrix-ness second,
    ## because a two-dimensional image is itself a matrix: without this order
    ## such an image would be misread as a raw transform to validate rather
    ## than an image whose transform is wanted
    if (isImage(x) || isImageGeometry(x) || !is.matrix(x))
    {
        geometry <- geometry(x)
        composeTransform(geometry@orientation, geometry@voxelSize)
    }
    else
        validateXform(x)
}

#' @rdname geometry
#' @export
`worldTransform<-` <- function (x, value)
{
    decomposed <- decomposeTransform(validateXform(value), spatial(x))
    updateGeometry(x, orientation = decomposed$orientation, voxelSize = decomposed$voxelSize)
}

#' @rdname geometry
#' @export
centre <- function (x)
{
    geometry <- geometry(x)
    n <- min(3L, length(geometry@dims))
    voxel <- c((geometry@dims[seq_len(n)] - 1) / 2, rep(0, 3L - n))
    as.vector(worldTransform(geometry) %*% c(voxel, 1))[1:3]
}

#' @rdname geometry
#' @export
center <- centre

#' @rdname geometry
#' @export
extent <- function (x)
{
    geometry <- geometry(x)
    geometry@dims * geometry@voxelSize
}

#' @rdname geometry
#' @export
sameGeometry <- function (x, y, tolerance = sqrt(.Machine$double.eps))
{
    x <- geometry(x)
    y <- geometry(y)
    if (!identical(x@dims, y@dims))
        return(FALSE)
    if (x@unit != "unknown" && y@unit != "unknown" && x@unit != y@unit)
        return(FALSE)
    isTRUE(all.equal(worldTransform(x), worldTransform(y), tolerance = tolerance, check.attributes = FALSE))
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
        if (is.null(bounds) && (isImage(x) || isImageGeometry(x) || !is.matrix(x)))
            bounds <- as.double(geometry(x)@dims[seq_len(min(3L, spatial(x)))])
        result <- roundPoints(result, round, bounds)
    }
    ## Standard one-based indexing for R
    if (!identical(type, "voxel"))
        result <- result + 1
    result[, seq_len(attr(points, "columns")), drop = FALSE]
}

#' @rdname geometry
#' @export
fromVoxel <- function (points, x, type = "world")
{
    points <- asPointMatrix(points)
    columns <- attr(points, "columns")
    ## The inverse of toVoxel()'s shift: points arrive one-based and are
    ## brought back to zero-based before reaching the affine, unless they are
    ## staying in voxel space, in which case there is nothing to shift. Only
    ## the columns supplied are shifted, so padding stays at the origin
    if (!identical(type, "voxel"))
        points[, seq_len(columns)] <- points[, seq_len(columns)] - 1
    result <- pointsFromVoxel(points, worldTransform(x), voxelSize(x), type)
    ## World space is always three-dimensional, even for a two-dimensional
    ## grid, which may be placed obliquely within it
    if (identical(type, "voxel"))
        result <- result[, seq_len(columns), drop = FALSE]
    result
}

## Points arrive as a vector (one point) or a matrix with one point per row.
## Points with fewer than three dimensions are padded with zeros, so that
## every dimensionality shares one code path, and the number supplied is remembered so that
## results can be trimmed back to match
asPointMatrix <- function (points)
{
    if (is.null(dim(points)))
        points <- matrix(points, nrow = 1L)
    points <- as.matrix(points)
    storage.mode(points) <- "double"
    columns <- ncol(points)
    if (columns < 1L || columns > 3L)
        stop("Points must be given as a matrix with one to three columns")
    if (columns < 3L)
        points <- cbind(points, matrix(0, nrow(points), 3L - columns))
    dimnames(points) <- NULL
    attr(points, "columns") <- columns
    points
}
