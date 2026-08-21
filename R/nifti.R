#' Interoperability with other image classes
#'
#' `imply` deliberately knows nothing about file formats: an image's geometry
#' is a 4x4 affine like any other. These bridges convert to and from the image
#' classes of packages that do handle formats, and are available whenever
#' those packages are installed.
#'
#' A NIfTI image's resolved xform maps voxel to world coordinates and bakes in
#' voxel size, which is exactly what `worldTransform()` returns for an imply
#' image, so the conversion is a copy of the data plus a transfer of geometry
#' (decomposed into `worldTransform()`'s rotation/translation and
#' `voxelSize()`'s scale on the way in).
#'
#' @param x An image to convert.
#' @param ... Further arguments to `denseImage()`.
#' @param datatype The NIfTI datatype to write, or `"auto"`.
#' @return `fromNifti()` and `fromMriImage()` return a `denseImage`.
#'   `toNifti()` returns an object of class `niftiImage`, and `toMriImage()`
#'   an `MriImage`.
#' @name interop
NULL

requirePackage <- function (name)
{
    if (!requireNamespace(name, quietly = TRUE))
        stop("The ", name, " package is needed for this conversion, but is not installed")
}

#' @rdname interop
#' @export
fromNifti <- function (x, ...)
{
    requirePackage("RNifti")

    values <- as.array(x)
    ## Strip the classes and bookkeeping attributes RNifti attaches, so what
    ## is left is a plain array
    for (name in c("class", ".nifti_image_ptr", ".nifti_image_ver", "pixdim", "pixunits", "imagedim"))
        attr(values, name) <- NULL

    dims <- dim(values)
    nSpatial <- min(3L, length(dims))
    units <- RNifti::pixunits(x)

    ## Voxel size and orientation both come from the resolved xform, rather
    ## than from the header's separately-stored pixdim field: RNifti::xform()
    ## already prefers the quaternion/qform-derived transform when available,
    ## which the NIfTI standard restricts to rigid+scale, and the two can
    ## genuinely disagree in files that have been edited by a tool that
    ## updates one and not the other
    image <- denseImage(values, spatial = nSpatial,
                        spaceUnit = if (length(units) > 0L) units[1L] else "unknown",
                        timeUnit = if (length(units) > 1L) units[2L] else "unknown",
                        ...)

    tryCatch(
        worldTransform(image) <- RNifti::xform(x)[1:4, 1:4],
        error = function (e) stop("NIfTI transform could not be interpreted as voxel geometry: ",
                                  conditionMessage(e), call. = FALSE))

    image
}

#' @rdname interop
#' @export
toNifti <- function (x, datatype = "auto")
{
    requirePackage("RNifti")

    if (isSparseImage(x) || isPackedImage(x))
        x <- asDense(x)
    x <- asDenseImage(x)

    ## RNifti exports no xform<-, and sform<- does not take on an image built
    ## from a bare array, so the transform goes in through the header fields it
    ## is actually stored in.
    ##
    ## Voxel dimensions have to travel in the same template rather than being
    ## assigned afterwards: RNifti's pixdim<- treats an assignment as a change
    ## of voxel size and rescales the transform to match, so setting it to the
    ## values worldTransform() already implies would double them
    transform <- worldTransform(x)
    voxelDims <- voxelSize(x)
    header <- rep(1, 8)
    header[seq_along(voxelDims) + 1L] <- voxelDims

    template <- list(pixdim = header,
                     sform_code = 2L,
                     srow_x = transform[1L, ],
                     srow_y = transform[2L, ],
                     srow_z = transform[3L, ])

    image <- RNifti::asNifti(as.array(x), reference = template, datatype = datatype)
    if (!identical(x@spaceUnit, "unknown"))
        RNifti::pixunits(image) <- c(x@spaceUnit, x@timeUnit)

    image
}

#' @rdname interop
#' @export
fromMriImage <- function (x, ...)
{
    requirePackage("tractor.base")

    values <- as.array(x)
    nSpatial <- min(3L, length(dim(values)))
    voxelDims <- abs(x$getVoxelDimensions())[seq_len(nSpatial)]

    transform <- x$getXform(implicit = TRUE)
    if (length(transform) == 0L || !identical(dim(transform), c(4L, 4L)))
        transform <- defaultXform(voxelDims)

    ## Orientation and voxel size both come from the transform, for the same
    ## reason as fromNifti(): a separately-tracked voxel size can disagree
    ## with the scale the transform implies
    denseImage(values, spatial = nSpatial, worldTransform = transform, ...)
}

#' @rdname interop
#' @export
toMriImage <- function (x)
{
    requirePackage("tractor.base")

    if (isSparseImage(x) || isPackedImage(x))
        x <- asDense(x)
    x <- asDenseImage(x)

    tractor.base::asMriImage(as.array(x), voxelDims = voxelSize(x), origin = originOf(x))
}

## An MriImage records the origin rather than the whole transform, so it is
## recovered as the voxel that maps to the world origin
originOf <- function (x)
{
    inverse <- invertXform(worldTransform(x))
    as.vector(inverse %*% c(0, 0, 0, 1))[1:3] + 1
}
