#' Image geometry
#'
#' Accessors for the geometry of an image: the number of spatial dimensions,
#' the voxel dimensions, the voxel-to-world transform and the anatomical
#' orientation it implies.
#'
#' These are deliberately free of any dependency on a file format. A NIfTI
#' image's transform is simply a 4x4 affine like any other.
#'
#' @param x An image, or for `orientation` either an image or a 4x4 matrix.
#' @param value A replacement value.
#' @param points A matrix of points, one per row and three columns.
#' @param type The coordinate convention of `points`: `"voxel"`, `"scaled"`
#'   (millimetres, ignoring rotation) or `"world"` (fully transformed).
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
spatial <- function (x) attr(x, "spatial") %||% min(3L, length(dim(x)))

#' @rdname geometry
#' @export
pixdim <- function (x) attr(x, "pixdim") %||% rep(1, spatial(x))

#' @rdname geometry
#' @export
`pixdim<-` <- function (x, value)
{
    value <- as.double(value)
    x <- asDenseImage(x)
    ## The transform encodes the voxel dimensions too, so it has to follow.
    ## Both assignments run the class validator
    x@pixdim <- value
    x@xform <- defaultXform(value)
    x
}

#' @rdname geometry
#' @export
xform <- function (x)
{
    if (is.matrix(x))
        validateXform(x)
    else
        attr(x, "xform") %||% defaultXform(pixdim(x))
}

#' @rdname geometry
#' @export
`xform<-` <- function (x, value)
{
    x <- asDenseImage(x)
    x@xform <- validateXform(value)
    x
}

defaultXform <- function (pixdim)
{
    result <- diag(4)
    n <- min(3L, length(pixdim))
    if (n > 0L)
        diag(result)[seq_len(n)] <- pixdim[seq_len(n)]
    result
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

#' @rdname geometry
#' @export
orientation <- function (x) orientationFromXform(xform(x))

#' @rdname geometry
#' @export
worldToVoxel <- function (points, x, type = "world", round = "none", bounds = NULL)
{
    points <- asPointMatrix(points)
    result <- pointsToVoxel(points, xform(x), pixdim(x), type)
    if (!identical(round, "none"))
    {
        if (is.null(bounds) && !is.matrix(x))
            bounds <- as.double(dim(x)[seq_len(min(3L, spatial(x)))])
        result <- roundPoints(result, round, bounds)
    }
    result
}

#' @rdname geometry
#' @export
voxelToWorld <- function (points, x, type = "world")
    pointsFromVoxel(asPointMatrix(points), xform(x), pixdim(x), type)

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
