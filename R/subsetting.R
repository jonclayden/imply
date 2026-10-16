## Indexing packed and sparse images. Both work from the same starting point:
## the subscripts are resolved to one-based linear indices in the image's own
## axis order, without materialising the image or an index array the size of
## it, and then mapped through the view to the storage. Extraction reads only
## the values asked for; replacement writes only them, into a copy, since R's
## copy-on-modify semantics apply as for any other object.
##
## A dense image is a plain array, and base R indexes it directly; see image.R

## The subscripts in a call, with a missing one standing for everything along
## its dimension, as TRUE does. Each is taken from `...` with ...elt(), which
## evaluates it in the caller's frame, as base R's own indexing would
captureSubscripts <- function (...)
{
    expressions <- as.list(substitute(list(...)))[-1L]
    lapply(seq_along(expressions), function (k) {
        if (identical(expressions[[k]], quote(expr = )))
            return(structure(TRUE, missing = TRUE))
        ## An image used as a subscript, such as x[x > 0], is its values
        i <- ...elt(k)
        if (isImage(i)) as.array(i) else i
    })
}

## One subscript along a dimension of extent `n`, as one-based indices. The
## rules are base R's, except that anything outside the dimension, and any
## missing index, is an error rather than giving or growing into NA
resolveSubscript <- function (i, n)
{
    if (is.null(i))
        return(integer(0))
    if (is.factor(i) || is.character(i))
        stop("Images can only be indexed by number or by logical value")
    if (is.logical(i))
    {
        if (length(i) > n)
            stop("Logical subscript is too long for a dimension of extent ", n)
        if (anyNA(i))
            stop("Missing values are not allowed in image subscripts")
        return(which(rep_len(i, n)))
    }
    if (!is.numeric(i))
        stop("Images can only be indexed by number or by logical value")
    if (anyNA(i))
        stop("Missing values are not allowed in image subscripts")

    i <- trunc(i)
    i <- i[i != 0]
    if (any(i < 0))
    {
        if (any(i > 0))
            stop("Positive and negative subscripts cannot be mixed")
        return(setdiff(seq_len(n), -i))
    }
    if (any(i > n))
        stop("Subscript out of bounds: ", max(i), " exceeds the extent ", n)
    i
}

## Linear indices, and the shape an extraction takes, for subscripts given in
## any of the forms base R accepts for an array: one per dimension, a single
## vector over all elements, or a matrix with one column per dimension
resolveIndices <- function (dims, subscripts)
{
    nDims <- length(dims)

    if (length(subscripts) == 1L)
    {
        i <- subscripts[[1L]]
        if (is.matrix(i) && is.numeric(i) && ncol(i) == nDims && nDims > 1L)
        {
            if (anyNA(i))
                stop("Missing values are not allowed in image subscripts")
            if (any(i < 1) || any(sweep(i, 2L, dims, ">")))
                stop("Subscript out of bounds")
            strides <- c(1, cumprod(as.double(dims)))[seq_len(nDims)]
            return(list(linear = as.vector((trunc(i) - 1) %*% strides) + 1, shape = NULL))
        }
        return(list(linear = as.double(resolveSubscript(i, prod(dims))), shape = NULL))
    }

    if (length(subscripts) != nDims)
        stop("The image has ", nDims, " dimensions, but ", length(subscripts), " subscripts were given")

    ## Each dimension's indices offset what came before, with the first
    ## varying fastest, as R orders an array
    perDimension <- mapply(resolveSubscript, subscripts, dims, SIMPLIFY = FALSE)
    stride <- 1
    linear <- 0
    for (k in seq_len(nDims))
    {
        linear <- as.vector(outer(linear, (perDimension[[k]] - 1) * stride, "+"))
        stride <- stride * dims[k]
    }
    list(linear = linear + 1, shape = lengths(perDimension))
}

## x[] is the whole array, dimensions and all, whatever drop says
isEverything <- function (subscripts)
    length(subscripts) == 1L && isTRUE(attr(subscripts[[1L]], "missing"))

## Extracted values take the shape of the subscripts, dropped as base R would
shapeExtracted <- function (values, shape, drop)
{
    if (is.null(shape))
        return(values)
    dim(values) <- shape
    if (drop) drop(values) else values
}

## The replacement values, recycled to the number of indices as base R
## recycles them
recycleValue <- function (value, n)
{
    if (isImage(value))
        value <- as.array(value)
    if (!is.atomic(value))
        stop("Replacement values must be atomic")
    value <- as.vector(value)
    if (n > 0L && length(value) == 0L)
        stop("Replacement has length zero")
    if (n > 0L && n %% length(value) != 0L)
        warning("The number of items to replace is not a multiple of the replacement length")
    rep_len(value, n)
}

## --- Packed images ---------------------------------------------------------

S7::method(`[`, packedImage) <- function (x, ..., drop = TRUE)
{
    subscripts <- captureSubscripts(...)
    if (isEverything(subscripts))
        return(as.array(x))
    indices <- resolveIndices(x@dims, subscripts)
    values <- narrowElements(x@values, x@storageType, prod(x@dims), storageIndices(x, indices$linear),
                             x@slope, x@intercept)
    shapeExtracted(values, indices$shape, drop)
}

## Values are stored under the image's existing type and scaling. Anything
## they cannot represent is refused, rather than the whole image being
## rescaled, which would change the rounding of every other value.
##
## Registered by calling method<- directly, rather than with `[<-`(x) <- value
## sugar: the sugar would create a namespace binding literally named `[<-`,
## which R CMD check's replacement function check chokes on. See image.R
S7::`method<-`(`[<-`, packedImage, function (x, ..., value)
{
    indices <- resolveIndices(x@dims, captureSubscripts(...))
    value <- recycleValue(value, length(indices$linear))
    if (is.complex(value))
        stop("Complex values cannot be stored in a packed image")

    x@values <- narrowAssign(x@values, x@storageType, prod(x@dims), storageIndices(x, indices$linear),
                             as.double(value), x@slope, x@intercept)
    x
})

## --- Sparse images ---------------------------------------------------------

S7::method(`[`, sparseImage) <- function (x, ..., drop = TRUE)
{
    subscripts <- captureSubscripts(...)
    if (isEverything(subscripts))
        return(as.array(x))
    indices <- resolveIndices(x@dims, subscripts)
    values <- sparseElements(x@mask, x@values, x@dims, spatial(x), storageIndices(x, indices$linear))
    shapeExtracted(values, indices$shape, drop)
}

## A non-zero value at a location not yet stored adds it; afterwards, any
## location left holding nothing but zeros is dropped, so the mask always
## says where the data are. Each assignment therefore makes a pass over the
## stored values, which is fine for occasional edits but not for setting many
## values one at a time in a loop: gather them and assign once. Values are
## promoted to the replacement's type where that is wider, as for an array
S7::`method<-`(`[<-`, sparseImage, function (x, ..., value)
{
    indices <- resolveIndices(x@dims, captureSubscripts(...))
    value <- recycleValue(value, length(indices$linear))
    if (length(value) == 0L)
        return(x)

    values <- x@values
    types <- c("logical", "integer", "double", "complex")
    wider <- types[max(match(c(typeof(values), typeof(value)), types))]
    if (is.na(wider))
        stop("Sparse images can only hold logical, integer, double or complex values")
    storage.mode(values) <- wider
    storage.mode(value) <- wider

    ## Stored order splits into a location and a position among the values
    ## held there, the spatial axes leading
    stored <- storageIndices(x, indices$linear) - 1
    nLocations <- locationCount(x)
    nElements <- elementCount(x)
    location <- stored %% nLocations + 1
    element <- stored %/% nLocations

    present <- maskToLogical(x@mask, nLocations)
    mask <- x@mask
    arriving <- !present[location] & (is.na(value) | value != 0)
    if (any(arriving))
    {
        present[location[arriving]] <- TRUE
        mask <- maskFromLogical(present)
        values <- repackValues(x@mask, values, mask, nLocations, nElements)
    }

    ## Zeros assigned where nothing is stored are already true
    rank <- cumsum(present)
    kept <- present[location]
    values[element[kept] + nElements * (rank[location[kept]] - 1) + 1] <- value[kept]

    rebuild(x, mask, values)
})
