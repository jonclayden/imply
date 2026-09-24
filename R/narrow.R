#' Narrow storage types
#'
#' R has no single-precision floating-point type, so a large image held as
#' `double` may cost twice the memory it needs, and potentially roughly twice
#' the time. Raw MRI is commonly 16-bit integer.
#'
#' Recognising this, a packed image stores its values using a narrow data type,
#' with optional affine scaling, so that an integer type can carry a range
#' it could not otherwise hold: `value = stored * slope + intercept`. When a
#' scaling is needed it is chosen automatically to map the data across the
#' whole of the type's range.
#'
#' The values live in a raw vector, so they are garbage-collected, serialise,
#' and survive a save and load like any other R object.
#'
#' @note Missing values can only be carried by `float32`; packing data
#' containing `NA` to an integer type is refused rather than silently losing
#' it. Note that `NA` and `NaN` are not distinguished once packed, since the
#' payload that separates them does not survive the narrowing.
#'
#' @param x An image or array.
#' @param type,storageType One of `"int8"`, `"uint8"`, `"int16"`, `"uint16"`,
#'   `"int32"` or `"float32"`.
#' @param slope,intercept Scaling applied to stored values. Chosen
#'   automatically when not given.
#' @param values A raw vector holding the packed values.
#' @param dims The full dimensions of the image.
#' @param spatial,voxelSize,worldTransform,unit,geometry Image geometry, as
#'   for [denseImage()].
#' @param layout How the image's axes map onto the order of `values`, as
#'   described for [storageLayout()]. The default is the identity, meaning
#'   that `values` are in the image's own axis order.
#' @param ... Further arguments to `denseImage()`.
#' @return An object of S7 class `packedImage` representing an image using a
#'   narrow, packed data representation, with properties corresponding to the
#'   arguments listed above.
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
        layout = S7::class_integer,
        geometry = imageGeometry
    ),
    validator = function (self) {
        if (length(self@storageType) != 1L || !self@storageType %in% storageTypes)
            return(paste0("@storageType must be one of ", paste(storageTypes, collapse = ", ")))
        if (length(self@slope) != 1L || is.na(self@slope) || self@slope == 0)
            return("@slope must be a single non-zero value")
        if (length(self@intercept) != 1L || is.na(self@intercept))
            return("@intercept must be a single value")

        if (anyNA(self@dims) || any(self@dims < 0L))
            return("@dims must not be missing or negative")
        mismatch <- geometryMismatch(self@geometry, self@dims) %||% checkLayout(self@layout, length(self@dims))
        if (!is.null(mismatch))
            return(mismatch)

        expected <- prod(self@dims) * storageTypeSize[[self@storageType]]
        if (length(self@values) != expected)
            return("@values is not the right length for the stated dimensions and storage type")

        NULL
    },
    constructor = function (values, storageType, dims, slope = 1, intercept = 0, spatial = NULL,
                            voxelSize = NULL, worldTransform = NULL, unit = NULL, geometry = NULL,
                            layout = NULL)
    {
        dims <- as.integer(dims)

        S7::new_object(S7::S7_object(),
            values = values,
            storageType = as.character(storageType),
            slope = as.double(slope),
            intercept = as.double(intercept),
            dims = dims,
            layout = as.integer(layout %||% seq_along(dims)),
            geometry = resolveGeometry(dims, spatial, voxelSize, worldTransform, unit, geometry))
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

    image <- asDense(x, ...)
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
                geometry = image@geometry)
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
    array(unpackNarrow(x@values, x@storageType, prod(x@dims), x@slope, x@intercept, x@dims, layoutArg(x)), x@dims)

S7::method(print, packedImage) <- function (x, ...)
{
    cat(sprintf("Packed image: %s (%s)\n", paste(x@dims, collapse = " x "), x@storageType))
    printGeometry(x@geometry)
    if (spatial(x) < length(x@dims))
        cat(sprintf("  Values per location: %d\n", prod(x@dims[-seq_len(spatial(x))])))
    if (x@slope != 1 || x@intercept != 0)
        cat(sprintf("  Scaling            : value = stored * %g + %g\n", x@slope, x@intercept))
    cat(sprintf("  Storage            : %s bytes, against %s as double\n",
                format(length(x@values), big.mark = ","),
                format(prod(x@dims) * 8, big.mark = ",")))
    printLayout(x@layout)

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
        return(narrowElements(x@values, x@storageType, prod(x@dims), storageIndices(x, i),
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
                denseImage(op(first, second), geometry = template@geometry)
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
