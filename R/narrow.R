#' Narrow storage types
#'
#' R has no single-precision type, so a large image held as `double` costs
#' twice the memory it needs, and — since these passes are bandwidth-bound —
#' roughly twice the time. Raw MRI is commonly 16-bit integer, four times
#' narrower again.
#'
#' A packed image stores its values in one of the NIfTI-style narrow types,
#' with an optional affine scaling, so that an integer type can carry a range
#' it could not otherwise hold: `value = stored * slope + intercept`. When a
#' scaling is needed it is chosen automatically to map the data across the
#' whole of the type's range.
#'
#' The values live in a raw vector, so they are garbage-collected, serialise,
#' and survive a save and load like any other R object. Nothing needs an
#' external pointer or a finalizer.
#'
#' Missing values can only be carried by `float32`; packing data containing
#' `NA` to an integer type is refused rather than silently losing it. Note
#' that `NA` and `NaN` are not distinguished once packed, since the payload
#' that separates them does not survive the narrowing.
#'
#' @param x An image or array.
#' @param type,storageType One of `"int8"`, `"uint8"`, `"int16"`, `"uint16"`,
#'   `"int32"` or `"float32"`.
#' @param slope,intercept Scaling applied to stored values. Chosen
#'   automatically when not given.
#' @param values A raw vector holding the packed values.
#' @param dims,spatial,pixdim,xform,spaceUnit,timeUnit Image geometry, as for
#'   [denseImage()].
#' @param template An image to take unspecified geometry from.
#' @param ... Further arguments to `denseImage()`.
#' @name packedImage
NULL

storageTypes <- c("int8", "uint8", "int16", "uint16", "int32", "float32")

storageTypeSize <- c(int8 = 1L, uint8 = 1L, int16 = 2L, uint16 = 2L, int32 = 4L, float32 = 4L)

#' @rdname packedImage
#' @export
packedImage <- S7::new_class("packedImage",
    properties = list(
        values = S7::class_raw,
        storageType = S7::class_character,
        slope = S7::class_double,
        intercept = S7::class_double,
        dims = S7::class_integer,
        spatial = S7::class_integer,
        pixdim = S7::class_double,
        xform = S7::class_double,
        spaceUnit = S7::class_character,
        timeUnit = S7::class_character
    ),
    validator = function (self) {
        nDims <- length(self@dims)

        if (length(self@storageType) != 1L || !self@storageType %in% storageTypes)
            return(paste0("@storageType must be one of ", paste(storageTypes, collapse = ", ")))
        if (length(self@slope) != 1L || is.na(self@slope) || self@slope == 0)
            return("@slope must be a single non-zero value")
        if (length(self@intercept) != 1L || is.na(self@intercept))
            return("@intercept must be a single value")

        if (length(self@spatial) != 1L || is.na(self@spatial))
            return("@spatial must be a single value")
        if (self@spatial < 0L || self@spatial > nDims)
            return(paste0("@spatial must be between 0 and ", nDims))

        expected <- prod(self@dims) * storageTypeSize[[self@storageType]]
        if (length(self@values) != expected)
            return("@values is not the right length for the stated dimensions and storage type")

        if (length(self@pixdim) != self@spatial)
            return(paste0("@pixdim must have one element per spatial dimension (", self@spatial, ")"))
        if (!identical(dim(self@xform), c(4L, 4L)))
            return("@xform must be a 4x4 matrix")

        NULL
    },
    constructor = function (values, storageType, dims, slope = 1, intercept = 0, spatial = NULL,
                            pixdim = NULL, xform = NULL, spaceUnit = NULL, timeUnit = NULL,
                            template = NULL)
    {
        dims <- as.integer(dims)
        nDims <- length(dims)
        spatial <- as.integer(spatial %||% attr(template, "spatial") %||% min(3L, nDims))

        pixdim <- as.double(pixdim %||% attr(template, "pixdim") %||% rep(1, max(spatial, 0L)))
        xform <- xform %||% attr(template, "xform") %||% defaultXform(pixdim)
        xform <- as.matrix(xform)
        storage.mode(xform) <- "double"
        dimnames(xform) <- NULL

        S7::new_object(S7::S7_object(),
            values = values,
            storageType = as.character(storageType),
            slope = as.double(slope),
            intercept = as.double(intercept),
            dims = dims,
            spatial = spatial,
            pixdim = pixdim,
            xform = xform,
            spaceUnit = as.character(spaceUnit %||% attr(template, "spaceUnit") %||% "unknown"),
            timeUnit = as.character(timeUnit %||% attr(template, "timeUnit") %||% "unknown"))
    })

S7::S4_register(packedImage)

#' @rdname packedImage
#' @export
isPackedImage <- function (x) S7::S7_inherits(x, packedImage)

#' @rdname packedImage
#' @export
asPacked <- function (x, type = "float32", slope = NULL, intercept = NULL, ...)
{
    if (isPackedImage(x) && identical(x@storageType, type))
        return(x)
    if (isSparseImage(x))
        x <- asDense(x)

    image <- asDenseImage(x, ...)
    values <- as.array(image)

    if (!typeof(values) %in% c("logical", "integer", "double"))
        stop("Only logical, integer and double data can be packed, not ", typeof(values))

    type <- match.arg(type, storageTypes)
    summary <- valueRange(values)

    if (summary$missing && type != "float32")
        stop("Data containing missing values can only be packed as float32, not ", type,
             ", which would lose them")

    if (is.null(slope) || is.null(intercept))
    {
        chosen <- calibrateStorage(type, summary$low, summary$high, summary$integral)
        slope <- slope %||% chosen$slope
        intercept <- intercept %||% chosen$intercept
    }

    packedImage(values = packNarrow(values, type, slope, intercept),
                storageType = type, dims = dim(image), slope = slope, intercept = intercept,
                spatial = image@spatial, pixdim = image@pixdim, xform = image@xform,
                spaceUnit = image@spaceUnit, timeUnit = image@timeUnit)
}

#' @rdname packedImage
#' @export
storageType <- function (x)
{
    if (isPackedImage(x))
        x@storageType
    else
        typeof(x)
}

S7::method(dim, packedImage) <- function (x) x@dims

S7::method(length, packedImage) <- function (x) prod(x@dims)

S7::method(as.array, packedImage) <- function (x, ...)
    array(unpackNarrow(x@values, x@storageType, prod(x@dims), x@slope, x@intercept), x@dims)

S7::method(print, packedImage) <- function (x, ...)
{
    cat(sprintf("Packed image: %s (%s)\n", paste(x@dims, collapse = " x "), x@storageType))
    if (x@spatial > 0L)
    {
        cat(sprintf("  Spatial dimensions : %s\n", paste(x@dims[seq_len(x@spatial)], collapse = " x ")))
        cat(sprintf("  Voxel dimensions   : %s %s\n",
                    paste(signif(x@pixdim, 4), collapse = " x "), x@spaceUnit))
        cat(sprintf("  Orientation        : %s\n", orientation(x)))
    }
    if (x@spatial < length(x@dims))
        cat(sprintf("  Values per location: %d\n", prod(x@dims[-seq_len(x@spatial)])))
    if (x@slope != 1 || x@intercept != 0)
        cat(sprintf("  Scaling            : value = stored * %g + %g\n", x@slope, x@intercept))
    cat(sprintf("  Storage            : %s bytes, against %s as double\n",
                format(length(x@values), big.mark = ","),
                format(prod(x@dims) * 8, big.mark = ",")))

    invisible(x)
}

## Indexing reads only the values asked for, rather than materialising the
## image. As elsewhere the result is a plain array
S7::method(`[`, packedImage) <- function (x, ..., drop = TRUE)
{
    indices <- as.list(substitute(list(...)))[-1L]
    supplied <- length(indices)
    absent <- vapply(indices, identical, NA, quote(expr = ))

    if (supplied == 0L || all(absent))
        return(as.array(x))

    if (supplied == 1L)
    {
        i <- ..1
        if (is.matrix(i) && ncol(i) == length(x@dims))
            i <- flattenIndices(array(0L, x@dims), i)
        else if (is.logical(i))
            i <- which(i)
        return(narrowElements(x@values, x@storageType, prod(x@dims), as.double(i),
                              x@slope, x@intercept))
    }

    call <- sys.call()
    call[[1L]] <- quote(`[`)
    call[[2L]] <- as.array(x)
    eval(call, parent.frame())
}

## Summaries are accumulated in double whatever the storage type, so the
## answer does not depend on how narrowly the values happen to be held
packedSummary <- function (x, na.rm = FALSE)
    narrowSummary(x@values, x@storageType, prod(x@dims), x@slope, x@intercept, na.rm)

registerPackedMethods <- function ()
{
    S7::`method<-`(base::sum, packedImage,
                   function (x, ..., na.rm = FALSE) packedSummary(x, na.rm)$sum)
    S7::`method<-`(base::min, packedImage,
                   function (x, ..., na.rm = FALSE) packedSummary(x, na.rm)$min)
    S7::`method<-`(base::max, packedImage,
                   function (x, ..., na.rm = FALSE) packedSummary(x, na.rm)$max)
    S7::`method<-`(base::range, packedImage,
                   function (x, ..., na.rm = FALSE) {
                       s <- packedSummary(x, na.rm)
                       c(s$min, s$max)
                   })
    S7::`method<-`(base::mean, packedImage,
                   function (x, ..., na.rm = FALSE) packedSummary(x, na.rm)$mean)

    ## Arithmetic materialises. Keeping a packed result would mean choosing a
    ## scaling for it, and the choice depends on the whole result, which has to
    ## be computed first anyway; pack the answer explicitly if it is wanted
    for (name in binaryOperators)
    {
        generic <- get(name, baseenv())
        handler <- local({
            op <- generic
            function (e1, e2) {
                first <- if (isPackedImage(e1)) as.array(e1) else asComparable(e1)
                second <- if (isPackedImage(e2)) as.array(e2) else asComparable(e2)
                template <- if (isPackedImage(e1)) e1 else e2
                denseImage(op(first, second), spatial = template@spatial, pixdim = template@pixdim,
                           xform = template@xform, spaceUnit = template@spaceUnit,
                           timeUnit = template@timeUnit)
            }
        })

        S7::`method<-`(generic, list(packedImage, packedImage), handler)
        S7::`method<-`(generic, list(packedImage, S7::class_numeric), handler)
        S7::`method<-`(generic, list(S7::class_numeric, packedImage), handler)
        S7::`method<-`(generic, list(packedImage, denseImage), handler)
        S7::`method<-`(generic, list(denseImage, packedImage), handler)
    }

    invisible(NULL)
}
