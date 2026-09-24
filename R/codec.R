#' Storage descriptors
#'
#' A storage descriptor says how the values of an image sit in a stream of
#' bytes: their type, byte order and position, any scaling applied to them, and
#' the order of the image's axes in the stream. It deliberately says nothing
#' about any particular file format. A format package parses a header into a
#' descriptor, and [readImageData()] and [writeImageData()] do the rest.
#'
#' Types are named by size and signedness: `"int8"`, `"uint8"`, `"int16"`,
#' `"uint16"`, `"int32"`, `"uint32"`, `"int64"`, `"float32"`, `"float64"`,
#' `"complex64"` and `"complex128"`, plus `"rgb24"` and `"rgba32"` for colour
#' and `"bit"` for one bit per value. The C-style names `"char"`, `"int"`,
#' `"float"` and `"double"` are accepted too, as `"int8"`, `"int32"`,
#' `"float32"` and `"float64"` respectively. C leaves the signedness of `char`
#' to the platform, so it is fixed as signed here; use `"uint8"` for unsigned
#' bytes.
#'
#' When read, integer types of up to 32 bits become R integers unless scaled,
#' and wider integers, floating-point types and scaled values become doubles.
#' 64-bit integers beyond 2^53 cannot be held exactly as doubles and are
#' refused. Colour values become integers from 0 to 255, with an extra
#' trailing dimension of three or four channels. Bits become logical values.
#' Scaling is not applied to complex or colour data.
#'
#' @param type The stored data type, as above.
#' @param endian The byte order of the stream: `"little"` or `"big"`.
#' @param offset The position of the first value, in bytes from the start of
#'   the stream.
#' @param slope,intercept Scaling, as `value = stored * slope + intercept`.
#'   When reading, `NULL` means no scaling. When writing, `NULL` means that
#'   a scaling is chosen to use the whole range of an integer type, if the
#'   data need one.
#' @param layout How the image's axes map onto the order of values in the
#'   stream, as described for [storageLayout()]. `NULL` means the stream
#'   holds the image in its own axis order, first axis fastest.
#' @param bitOrder For `"bit"` data, whether the first value of each byte is
#'   its least or most significant bit.
#' @return An object of S7 class `storageDescriptor`.
#' @examples
#' storageDescriptor("int16", endian = "big", offset = 352)
#' @export
storageDescriptor <- S7::new_class("storageDescriptor",
    properties = list(
        type = S7::class_character,
        endian = S7::class_character,
        offset = S7::class_double,
        slope = S7::class_double,
        intercept = S7::class_double,
        layout = S7::class_integer,
        bitOrder = S7::class_character
    ),
    validator = function (self) {
        if (length(self@type) != 1L || !self@type %in% names(codecTypeSizes))
            return(paste0("@type must be one of ", paste(names(codecTypeSizes), collapse = ", ")))
        if (length(self@endian) != 1L || !self@endian %in% c("little", "big"))
            return("@endian must be \"little\" or \"big\"")
        if (length(self@offset) != 1L || is.na(self@offset) || self@offset < 0)
            return("@offset must be a single non-negative value")
        if (length(self@slope) != 1L || length(self@intercept) != 1L)
            return("@slope and @intercept must each be a single value, which may be NA")
        if (isTRUE(self@slope == 0) || isTRUE(!is.finite(self@slope) && !is.na(self@slope)))
            return("@slope must be finite and non-zero")
        if (length(self@layout) > 0L)
        {
            problem <- checkLayout(self@layout, length(self@layout))
            if (!is.null(problem))
                return(problem)
        }
        if (length(self@bitOrder) != 1L || !self@bitOrder %in% c("lsb", "msb"))
            return("@bitOrder must be \"lsb\" or \"msb\"")
        NULL
    },
    constructor = function (type, endian = "little", offset = 0, slope = NULL, intercept = NULL,
                            layout = NULL, bitOrder = "lsb")
    {
        type <- as.character(type)
        if (length(type) == 1L && type %in% names(codecTypeAliases))
            type <- codecTypeAliases[[type]]
        S7::new_object(S7::S7_object(),
            type = type,
            endian = as.character(endian),
            offset = as.double(offset),
            slope = as.double(slope %||% NA_real_),
            intercept = as.double(intercept %||% NA_real_),
            layout = as.integer(layout %||% integer(0)),
            bitOrder = as.character(bitOrder))
    })

S7::S4_register(storageDescriptor)

S7::method(print, storageDescriptor) <- function (x, ...)
{
    cat(sprintf("Storage descriptor: %s, %s-endian, from byte %s\n", x@type, x@endian, format(x@offset)))
    if (!is.na(x@slope) || !is.na(x@intercept))
        cat(sprintf("  Scaling: value = stored * %g + %g\n", x@slope %|NA|% 1, x@intercept %|NA|% 0))
    if (length(x@layout) > 0L && !isIdentityLayout(x@layout))
        cat(sprintf("  Layout : %s\n", paste(x@layout, collapse = ", ")))
    invisible(x)
}

`%|NA|%` <- function (x, y) if (is.na(x)) y else x

codecTypeSizes <- c(int8 = 1, uint8 = 1, int16 = 2, uint16 = 2, int32 = 4, uint32 = 4, int64 = 8,
                    float32 = 4, float64 = 8, complex64 = 8, complex128 = 16, rgb24 = 3, rgba32 = 4,
                    bit = 1/8)

codecTypeAliases <- c(char = "int8", int = "int32", float = "float32", double = "float64")

## Colour is stored as bytes with the channel varying fastest, which is an
## extra storage axis of rank one in front of the image's own. The image sees
## the channels as a trailing value dimension
colourChannels <- c(rgb24 = 3L, rgba32 = 4L)

#' Reading and writing image data
#'
#' `readImageData()` reads the values of an image from a connection, as
#' described by a [storageDescriptor()], and returns an image. `writeImageData()`
#' writes them. Neither knows anything about file formats, which are the
#' business of packages built on this one: a format package parses its header,
#' builds a descriptor and a geometry, and calls these to handle the data.
#'
#' Bytes are read through R connections, so compressed streams need nothing
#' special: [gzfile()] reads gzipped and plain files alike. They are read in
#' chunks, so the only large allocation is the result itself. A gzipped stream
#' cannot be read backwards, so volumes are always read in the order they are
#' stored, whatever order they are asked for in.
#'
#' The image comes back dense, packed or sparse, as `as` says. The default,
#' `"auto"`, gives a sparse image if a `mask` is given; otherwise a packed image
#' if the stored type is one a [packedImage] can hold, keeping the file's own
#' narrow representation and scaling; and otherwise a dense image. The class
#' returned therefore depends on the data, so code that relies on one should
#' ask for it.
#'
#' A packed image keeps the stream's own order, with the layout recording how
#' its axes map onto it, so no data is moved; see [reorient()]. A dense image
#' is always in its own axis order, so a stream in any other order is
#' rearranged as it is read.
#'
#' A sparse image keeps exactly the locations selected by `mask`, if one is
#' given, and otherwise those where any value is non-zero. Only volume-ordered
#' data can be read straight into sparse form, without the dense array ever
#' existing; other layouts, and bit data, are read densely first.
#'
#' @param con A connection, or the path to a file. A connection must be open
#'   for binary reading or writing. If its position can be queried, it is
#'   moved forwards to the descriptor's offset, or when writing, padded with
#'   zeros up to it; otherwise it is taken to be already there. A path is
#'   opened, read or written from the start, and closed.
#' @param descriptor A [storageDescriptor()].
#' @param dims The full dimensions of the image, in its own axis order.
#' @param geometry An [imageGeometry][geometry], or an image to take one from,
#'   for the spatial grid. By default the leading three dimensions (or all of
#'   them, if fewer) are spatial, with unit voxels at the origin.
#' @param as The representation wanted: `"auto"`, `"dense"`, `"sparse"` or
#'   `"packed"`.
#' @param volumes Optionally, which volumes to read, as one-based indices into
#'   the values held at each location, taken together. The result then has the
#'   spatial dimensions followed by one dimension over the volumes read, in the
#'   order given. Only data stored a volume at a time can be read this way.
#' @param mask Optionally, the spatial locations to keep, in any form
#'   [asMaskVector()] accepts. Only meaningful for a sparse result.
#' @param image An image or array to write.
#' @return `readImageData()` returns an image. `writeImageData()` returns,
#'   invisibly, the descriptor as used, with any scaling it chose filled in,
#'   which is what a format package needs to write its header.
#' @examples
#' image <- denseImage(array(1:24, c(2, 3, 4)))
#' descriptor <- storageDescriptor("int16", endian = "big")
#' file <- tempfile()
#' writeImageData(image, file, descriptor)
#' readImageData(file, descriptor, dim(image))
#' @export
readImageData <- function (con, descriptor, dims, geometry = NULL, as = c("auto", "dense", "sparse", "packed"),
                           volumes = NULL, mask = NULL)
{
    if (!S7::S7_inherits(descriptor, storageDescriptor))
        stop("A storageDescriptor is needed to describe the data")
    as <- match.arg(as)
    dims <- as.integer(dims)
    geometry <- resolveGeometry(dims, from = geometry)
    nSpatial <- length(geometry@dims)
    type <- descriptor@type
    slope <- descriptor@slope %|NA|% 1
    intercept <- descriptor@intercept %|NA|% 0
    layout <- if (length(descriptor@layout) == 0L) seq_along(dims) else descriptor@layout
    if (length(layout) != length(dims))
        stop("The descriptor's layout has ", length(layout), " entries, but the image has ", length(dims), " dimensions")

    ## Colour becomes bytes with an extra, fastest-varying storage axis, seen by
    ## the image as a trailing dimension
    streamDims <- dims
    if (type %in% names(colourChannels))
    {
        streamDims <- c(dims, colourChannels[[type]])
        layout <- c(sign(layout) * (abs(layout) + 1L), 1L)
        type <- "uint8"
        slope <- 1
        intercept <- 0
    }

    if (as == "auto")
        as <- if (!is.null(mask)) "sparse" else if (descriptor@type %in% storageTypes) "packed" else "dense"
    if (!is.null(mask) && as != "sparse")
        stop("A mask is only used when reading to a sparse image")
    if (as == "packed" && !type %in% storageTypes)
        stop("Data of type ", descriptor@type, " cannot be kept packed; read them densely instead")

    ## Whether values are stored a volume at a time: the spatial axes first, in
    ## any order and direction among themselves, then the rest in order
    blockOrdered <- all(abs(layout[seq_len(nSpatial)]) <= nSpatial) &&
                    isIdentityLayout(layout[-seq_len(nSpatial)] - nSpatial)
    nVolumes <- prod(streamDims[-seq_len(nSpatial)])

    blocks <- NULL
    outDims <- streamDims
    if (!is.null(volumes))
    {
        volumes <- as.integer(volumes)
        if (length(volumes) == 0L || anyNA(volumes) || any(volumes < 1L) || any(volumes > nVolumes))
            stop("Volumes must be between 1 and ", nVolumes)
        if (!blockOrdered)
            stop("Volumes can only be selected from data stored a volume at a time")
        blocks <- as.double(volumes - 1L)
        outDims <- c(streamDims[seq_len(nSpatial)], length(volumes))
    }

    ## Sparse results are built a volume at a time, which needs volume order.
    ## Anything else is read densely and converted
    viaDense <- (as == "sparse" && (!blockOrdered || type == "bit"))
    selected <- if (is.null(mask)) NULL else asMaskVector(mask, geometry@dims)

    if (is.character(con))
    {
        con <- gzfile(con, "rb")
        on.exit(close(con))
    }
    reader <- streamReader(con)
    reader$skip(descriptor@offset - reader$position())

    ## Reads that keep the stream's order carry its layout with them, reduced
    ## to the spatial part plus one volume axis when volumes were selected
    keptLayout <- if (is.null(blocks)) layout else c(layout[seq_len(nSpatial)], nSpatial + 1L)
    if (as == "sparse" && !viaDense && is.null(blocks))
        keptLayout <- c(layout[seq_len(nSpatial)], seq_len(length(outDims) - nSpatial) + nSpatial)

    result <- decodeImage(reader$read, reader$skip, type, descriptor@endian != .Platform$endian,
                          slope, intercept, streamDims, if (isIdentityLayout(layout)) NULL else layout,
                          nSpatial, if (viaDense) "dense" else as, blocks,
                          if (viaDense) NULL else selected, descriptor@bitOrder == "msb")

    if (as == "packed")
        return(packedImage(values = result, storageType = type, dims = outDims, slope = slope,
                           intercept = intercept, geometry = geometry, layout = keptLayout))

    if (as == "sparse" && !viaDense)
    {
        ## Sparse storage has one row per location, one column per value held
        values <- result$values
        values <- if (length(outDims) > nSpatial) values else as.vector(values)
        return(sparseImage(mask = result$mask, values = values, dim = outDims, geometry = geometry,
                           layout = keptLayout))
    }

    dim(result) <- outDims
    image <- denseImage(result, geometry = geometry)
    if (as == "sparse")
        image <- sparseFromDense(image, selected)
    image
}

## A sparse image from a dense one, keeping exactly the locations selected if
## a selection is given, and otherwise the non-zero ones
sparseFromDense <- function (image, selected = NULL)
{
    if (is.null(selected))
        return(asSparse(image))

    nSpatial <- spatial(image)
    values <- as.array(image)
    nLocations <- prod(dim(image)[seq_len(nSpatial)])
    dim(values) <- c(nLocations, length(values) %/% max(nLocations, 1L))
    sparseImage(mask = selected, values = t(values[selected, , drop = FALSE]), dim = dim(image),
                geometry = image@geometry)
}

## Closures over a connection for the compiled decoder. Plain files are
## skipped through by seeking; anything else, including compressed streams,
## is read and discarded, which is the only way forward through a gzip stream
streamReader <- function (con)
{
    canSeek <- isTRUE(tryCatch(isSeekable(con), error = function (e) FALSE)) &&
               identical(summary(con)$class, "file")

    read <- function (n) readBin(con, "raw", n)

    skip <- function (n)
    {
        if (n < 0)
            stop("The connection is already past the start of the data")
        if (n == 0)
            return(invisible(NULL))
        if (canSeek)
            seek(con, n, origin = "current")
        else
        {
            while (n > 0)
            {
                chunk <- min(n, 2^26)
                if (length(readBin(con, "raw", chunk)) < chunk)
                    stop("Unexpected end of data while skipping forward")
                n <- n - chunk
            }
        }
        invisible(NULL)
    }

    position <- function ()
    {
        where <- tryCatch(seek(con), error = function (e) NA_real_)
        if (is.na(where)) 0 else where
    }

    list(read = read, skip = skip, position = position)
}

#' @rdname readImageData
#' @export
writeImageData <- function (image, con, descriptor)
{
    if (!S7::S7_inherits(descriptor, storageDescriptor))
        stop("A storageDescriptor is needed to describe the data")

    type <- descriptor@type
    dims <- dim(image) %||% length(image)
    swap <- (descriptor@endian != .Platform$endian)

    if (type %in% names(colourChannels))
    {
        ## The channels are the image's last dimension, but the stream's
        ## fastest-varying axis; the descriptor's layout covers the rest
        channels <- colourChannels[[type]]
        if (utils::tail(dims, 1L) != channels)
            stop("Colour data of type ", type, " need a last dimension of ", channels)
        layout <- if (length(descriptor@layout) == 0L) seq_len(length(dims) - 1L) else descriptor@layout
        if (length(layout) != length(dims) - 1L)
            stop("The descriptor's layout has ", length(layout), " entries, but the image has ",
                 length(dims) - 1L, " dimensions besides its colour channels")
        layout <- c(sign(layout) * (abs(layout) + 1L), 1L)
        type <- "uint8"
        descriptor@slope <- descriptor@intercept <- NA_real_
    }
    else
        layout <- if (length(descriptor@layout) == 0L) seq_along(dims) else descriptor@layout

    if (length(layout) != length(dims))
        stop("The descriptor's layout has ", length(layout), " entries, but the image has ", length(dims), " dimensions")

    ## Packed data already in the form wanted are written as they stand
    if (isPackedImage(image) && image@storageType == type &&
        (is.na(descriptor@slope) || isTRUE(all.equal(descriptor@slope, image@slope))) &&
        (is.na(descriptor@intercept) || isTRUE(all.equal(descriptor@intercept, image@intercept))) &&
        identical(as.integer(layout), image@layout))
    {
        descriptor@slope <- image@slope
        descriptor@intercept <- image@intercept
        bytes <- if (swap) swapBytes(image@values, storageTypeSize[[type]]) else image@values
    }
    else
    {
        values <- as.array(image)
        if (is.na(descriptor@slope) || is.na(descriptor@intercept))
        {
            chosen <- chooseScaling(values, type)
            descriptor@slope <- descriptor@slope %|NA|% chosen$slope
            descriptor@intercept <- descriptor@intercept %|NA|% chosen$intercept
        }
        bytes <- encodeImage(values, type, swap, descriptor@slope, descriptor@intercept, dims,
                             if (isIdentityLayout(layout)) NULL else layout, descriptor@bitOrder == "msb")
    }

    if (is.character(con))
    {
        con <- if (grepl("\\.gz$", con, ignore.case = TRUE)) gzfile(con, "wb") else file(con, "wb")
        on.exit(close(con))
    }
    position <- tryCatch(seek(con), error = function (e) NA_real_)
    if (!is.na(position))
    {
        if (position > descriptor@offset)
            stop("The connection is already past the offset the data should start at")
        if (position < descriptor@offset)
            writeBin(raw(descriptor@offset - position), con)
    }
    writeBin(bytes, con)

    if (descriptor@type %in% names(colourChannels))
        descriptor@slope <- descriptor@intercept <- NA_real_
    invisible(descriptor)
}

## A scaling for writing: the narrow integer types get one that uses their
## whole range when the data need it, as for asPacked(); everything else is
## written as it stands
chooseScaling <- function (values, type)
{
    if (!type %in% setdiff(storageTypes, "float32") || !typeof(values) %in% c("logical", "integer", "double"))
        return(list(slope = 1, intercept = 0))
    summary <- valueRange(values)
    if (!is.finite(summary$low) || !is.finite(summary$high))
        return(list(slope = 1, intercept = 0))
    calibrateStorage(type, summary$low, summary$high, summary$integral)
}
