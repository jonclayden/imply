## The storage codec: reading and writing image data described by a storage
## descriptor, with no file format involved. Everything is a round trip,
## checked against the array that went in, across types, byte orders, layouts
## and representations.

set.seed(14)
dims <- c(4L, 5L, 3L, 2L)
signed <- array(sample(-50:50, prod(dims), replace = TRUE), dims)
signed[abs(signed) < 20] <- 0L
unsigned <- abs(signed)

integerTypes <- c("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64")
floatTypes <- c("float32", "float64")
packable <- c("int8", "uint8", "int16", "uint16", "int32", "float32")

## Layouts: the identity, a volume-ordered one with a spatial axis reversed
## and two swapped, and an interleaved one with the volume axis fastest
layouts <- list(identity = NULL, reordered = c(2L, -1L, 3L, 4L), interleaved = c(2L, 3L, 4L, 1L))

## --- Round trips -----------------------------------------------------------

for (type in c(integerTypes, floatTypes)) for (endian in c("little", "big")) for (name in names(layouts))
{
    info <- paste(type, endian, name)
    source <- if (startsWith(type, "u")) unsigned else signed
    file <- tempfile(fileext = if (endian == "big") ".gz" else "")
    descriptor <- storageDescriptor(type, endian = endian, offset = 17, layout = layouts[[name]])
    used <- writeImageData(source, file, descriptor)

    expect_equal(used@slope, 1, info = info)
    for (as in c("dense", "sparse", if (type %in% packable) "packed"))
    {
        image <- readImageData(file, used, dims, as = as)
        expect_equal(as.vector(as.array(image)), as.vector(source) * 1, info = paste(info, as))
        expect_equal(dim(image), dims, info = paste(info, as))
    }

    ## Unscaled integer types up to 32 bits come back as integers
    dense <- readImageData(file, used, dims, as = "dense")
    expect_identical(typeof(dense), if (type %in% c("uint32", "int64") || type %in% floatTypes) "double" else "integer",
                     info = info)

    if (name != "interleaved")
    {
        ## Selected volumes, in any order and with repeats
        subset <- readImageData(file, used, dims, as = "dense", volumes = c(2, 1, 2))
        expect_equal(as.vector(as.array(subset)), as.vector(source[, , , c(2, 1, 2)]) * 1, info = info)
        expect_equal(dim(subset), c(4L, 5L, 3L, 3L), info = info)
        if (type %in% packable)
            expect_equal(as.vector(as.array(readImageData(file, used, dims, as = "packed", volumes = 2))),
                         as.vector(source[, , , 2]) * 1, info = info)

        ## A mask gives a sparse image holding exactly its locations
        brain <- array(runif(60) < 0.5, dims[1:3])
        sparse <- readImageData(file, used, dims, mask = brain)
        expect_true(isSparseImage(sparse), info = info)
        expect_identical(mask(sparse), brain, info = info)
        expect_equal(as.vector(as.array(sparse)), as.vector(source * as.vector(brain)) * 1, info = info)
    }
    else
        expect_error(readImageData(file, used, dims, volumes = 1), "a volume at a time", info = info)
}

## A stream in another order is a view when kept packed, with no data moved,
## and is rearranged when read densely
file <- tempfile()
descriptor <- storageDescriptor("int16", layout = c(2L, -1L, 3L, 4L))
writeImageData(signed, file, descriptor)
packed <- readImageData(file, descriptor, dims, as = "packed")
expect_identical(storageLayout(packed), c(2L, -1L, 3L, 4L))
expect_identical(packed@values, readBin(file, "raw", 2 * prod(dims)))
expect_identical(as.array(readImageData(file, descriptor, dims, as = "dense")), signed)

## The file really is in the order described: the first value stored is the
## last along the reversed first axis
expect_equal(readBin(file, "integer", 1L, size = 2L), signed[4, 1, 1, 1])

## --- The automatic choice of representation --------------------------------

file <- tempfile()
writeImageData(signed, file, storageDescriptor("int16"))
expect_true(isPackedImage(readImageData(file, storageDescriptor("int16"), dims)))
expect_true(isSparseImage(readImageData(file, storageDescriptor("int16"), dims, mask = signed[, , , 1] != 0)))
writeImageData(signed, file, storageDescriptor("float64"))
expect_true(isDenseImage(readImageData(file, storageDescriptor("float64"), dims)))
expect_error(readImageData(file, storageDescriptor("float64"), dims, as = "packed"), "cannot be kept packed")
expect_error(readImageData(file, storageDescriptor("float64"), dims, as = "dense", mask = signed[, , , 1] != 0),
             "only used when reading to a sparse")

## --- Scaling ---------------------------------------------------------------

smooth <- array(seq(-1, 1, length.out = 60), c(3L, 4L, 5L))
file <- tempfile()

## Chosen automatically on writing, using the type's whole range, and returned
used <- writeImageData(smooth, file, storageDescriptor("int16"))
expect_true(used@slope != 1)
expect_true(max(abs(as.array(readImageData(file, used, dim(smooth))) - smooth)) < 1e-4)
expect_true(isPackedImage(readImageData(file, used, dim(smooth))))
expect_equal(readImageData(file, used, dim(smooth))@slope, used@slope)

## ...or given, in which case values outside the type's range are refused
expect_error(writeImageData(smooth * 1e6, file, storageDescriptor("int16", slope = 1, intercept = 0)), "outside the range")

## Scaled integers are read as doubles
expect_identical(typeof(as.array(readImageData(file, used, dim(smooth), as = "dense"))), "double")

## Missing values need a floating-point type
withNA <- replace(smooth, 5, NA)
expect_error(writeImageData(withNA, file, storageDescriptor("int16")), "Missing values")
writeImageData(withNA, file, storageDescriptor("float32"))
expect_true(is.na(as.array(readImageData(file, storageDescriptor("float32"), dim(smooth)))[5]))
writeImageData(withNA, file, storageDescriptor("float64"))
expect_identical(as.array(readImageData(file, storageDescriptor("float64"), dim(smooth))), withNA)

## Packed data already in the form wanted are written byte for byte
packed <- asPacked(denseImage(smooth), "int16")
used <- writeImageData(packed, file, storageDescriptor("int16", endian = "big"))
back <- readImageData(file, used, dim(smooth))
expect_identical(back@values, packed@values)
expect_equal(back@slope, packed@slope)

## --- Other types -----------------------------------------------------------

## Colour: channels are interleaved in the stream, and a trailing dimension in
## the image
colour <- array(sample(0:255, 4 * 5 * 3, replace = TRUE), c(4L, 5L, 3L))
file <- tempfile()
writeImageData(colour, file, storageDescriptor("rgb24"))
expect_identical(readBin(file, "raw", 3L), as.raw(colour[1, 1, ]))
expect_identical(as.array(readImageData(file, storageDescriptor("rgb24"), c(4L, 5L))), colour)
expect_error(writeImageData(colour[, , 1:2], file, storageDescriptor("rgb24")), "last dimension of 3")
withAlpha <- array(sample(0:255, 2 * 2 * 4, replace = TRUE), c(2L, 2L, 4L))
writeImageData(withAlpha, file, storageDescriptor("rgba32"))
expect_identical(as.array(readImageData(file, storageDescriptor("rgba32"), c(2L, 2L))), withAlpha)

## Bits, in either order within a byte, including volumes that start part way
## through one
bits <- array(runif(3 * 3 * 1 * 5) < 0.5, c(3L, 3L, 1L, 5L))
for (order in c("lsb", "msb"))
{
    descriptor <- storageDescriptor("bit", bitOrder = order)
    writeImageData(bits, file, descriptor)
    expect_equal(file.size(file), ceiling(length(bits) / 8))
    expect_identical(as.array(readImageData(file, descriptor, dim(bits))), bits, info = order)
    expect_identical(as.vector(as.array(readImageData(file, descriptor, dim(bits), volumes = c(4, 2)))),
                     as.vector(bits[, , , c(4, 2)]), info = order)
}
expect_true(isSparseImage(readImageData(file, storageDescriptor("bit", bitOrder = "msb"), dim(bits), as = "sparse")))

## The order really is different on disk: the first value is the lowest bit of
## the first byte, or the highest
first <- c(TRUE, rep(FALSE, 7L))
writeImageData(first, file, storageDescriptor("bit", bitOrder = "lsb"))
expect_identical(readBin(file, "raw", 1L), as.raw(0x01))
writeImageData(first, file, storageDescriptor("bit", bitOrder = "msb"))
expect_identical(readBin(file, "raw", 1L), as.raw(0x80))

## Complex, in both widths
phase <- array(complex(real = rnorm(12), imaginary = rnorm(12)), c(3L, 4L))
writeImageData(phase, file, storageDescriptor("complex128", endian = "big"))
expect_identical(as.array(readImageData(file, storageDescriptor("complex128", endian = "big"), dim(phase))), phase)
writeImageData(phase, file, storageDescriptor("complex64"))
expect_equal(as.array(readImageData(file, storageDescriptor("complex64"), dim(phase))), phase, tolerance = 1e-6)

## 64-bit integers beyond 2^53 cannot be doubles exactly, so are refused
connection <- file(file, "wb")
writeBin(c(0L, 4194304L), connection, size = 4L)
close(connection)
expect_error(readImageData(file, storageDescriptor("int64"), 1L), "2\\^53")

## R's NA integer cannot come from an int32 file; packed storage can hold it
connection <- file(file, "wb")
writeBin(c(1L, NA_integer_), connection, size = 4L)
close(connection)
expect_error(readImageData(file, storageDescriptor("int32"), 2L, as = "dense"), "packed instead")
expect_equal(as.vector(as.array(readImageData(file, storageDescriptor("int32"), 2L, as = "packed"))), c(1, -2147483648))

## --- Descriptors -----------------------------------------------------------

## C-style names, at fixed sizes
expect_equal(storageDescriptor("char")@type, "int8")
expect_equal(storageDescriptor("int")@type, "int32")
expect_equal(storageDescriptor("float")@type, "float32")
expect_equal(storageDescriptor("double")@type, "float64")

expect_error(storageDescriptor("int12"), "@type must be one of")
expect_error(storageDescriptor("int16", endian = "middle"), "@endian")
expect_error(storageDescriptor("int16", offset = -1), "@offset")
expect_error(storageDescriptor("int16", slope = 0), "@slope")
expect_error(storageDescriptor("int16", layout = c(1L, 1L)), "exactly once")
expect_error(readImageData(file, list(type = "int16"), 2L), "storageDescriptor is needed")
expect_error(readImageData(file, storageDescriptor("int16", layout = 1:3), 2L), "layout has 3 entries")
expect_true(any(grepl("big-endian", capture.output(print(storageDescriptor("int16", endian = "big"))))))

## --- Connections and positions ---------------------------------------------

## An open connection is moved forward to the offset, so a header can be read
## first and the same connection handed over
connection <- file(file, "wb")
writeBin(charToRaw("HEADER"), connection)
writeImageData(signed, connection, storageDescriptor("int32", offset = 10))
close(connection)
expect_equal(file.size(file), 10 + 4 * length(signed))

connection <- file(file, "rb")
invisible(readBin(connection, "raw", 6L))
image <- readImageData(connection, storageDescriptor("int32", offset = 10), dims, as = "dense")
close(connection)
expect_identical(as.array(image), signed)

## The same through a gzip stream, which can only be read forwards
file <- tempfile(fileext = ".gz")
connection <- gzfile(file, "wb")
writeBin(charToRaw("HEADER"), connection)
writeImageData(signed, connection, storageDescriptor("int32", offset = 10))
close(connection)
connection <- gzfile(file, "rb")
invisible(readBin(connection, "raw", 6L))
image <- readImageData(connection, storageDescriptor("int32", offset = 10), dims, as = "dense", volumes = 2)
close(connection)
expect_identical(as.vector(as.array(image)), as.vector(signed[, , , 2]))

## A stream that ends early is reported rather than padded
expect_error(readImageData(file, storageDescriptor("int32", offset = 10), c(dims, 2L)), "Unexpected end of data")

## The geometry given is used for the image, and must match its grid
grid <- imageGeometry(c(4L, 5L, 3L), voxelSize = c(2, 2, 3), unit = "mm")
image <- readImageData(file, storageDescriptor("int32", offset = 10), dims, geometry = grid)
expect_identical(geometry(image), grid)
expect_error(readImageData(file, storageDescriptor("int32", offset = 10), dims,
                           geometry = imageGeometry(c(4L, 5L, 2L))), "grid of 4 x 5 x 2")
