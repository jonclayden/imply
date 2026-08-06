## Arithmetic, comparison and summary on sparse images.
##
## The governing question for every operation is what it does to an absent
## location, whose values are all implicitly zero. If the answer is still zero
## the operation can be applied to the stored values alone and the result stays
## sparse. If it is not -- as for `x + 1`, where every absent location becomes
## one -- then the result genuinely is dense, and saying so is more honest than
## re-sparsifying at some threshold and hoping. tractor.base takes the latter
## course, densifying for all arithmetic and re-packing only when the result is
## at least 75% zeros.

binaryOperators <- c("+", "-", "*", "/", "^", "%%", "%/%",
                     "==", "!=", "<", ">", "<=", ">=", "&", "|")

## Functions that leave zero alone are safe to push through the stored values;
## the rest have to densify first
mathFunctions <- c("sqrt", "abs", "floor", "ceiling", "trunc", "round", "signif",
                   "sin", "tan", "sinh", "tanh", "asin", "atan", "expm1", "log1p", "sign",
                   "exp", "log", "log2", "log10", "cos", "cosh", "acos", "gamma", "lgamma")

summaryFunctions <- c("sum", "max", "min", "range", "any", "all", "prod")

zeroOf <- function (x) vector(typeof(x), 1L)

## A length-one probe answering "what does this operation make of a location
## that holds nothing?"
absentResult <- function (op, first, second)
    tryCatch(op(first, second), error = function (e) NULL)

isStillAbsent <- function (value)
{
    !is.null(value) && length(value) == 1L && !is.na(value) &&
        (isFALSE(value) || isTRUE(all.equal(as.vector(value), 0)))
}

geometryOf <- function (x)
    list(spatial = x@spatial, pixdim = x@pixdim, xform = x@xform,
         spaceUnit = x@spaceUnit, timeUnit = x@timeUnit)

rebuild <- function (template, mask, values)
{
    tight <- tightenMask(mask, values, locationCount(template), elementCount(template))
    sparseImage(mask = tight$mask, values = tight$values, dim = template@dims,
                spatial = template@spatial, pixdim = template@pixdim, xform = template@xform,
                spaceUnit = template@spaceUnit, timeUnit = template@timeUnit)
}

denseFrom <- function (template, values)
    denseImage(array(values, template@dims), spatial = template@spatial, pixdim = template@pixdim,
               xform = template@xform, spaceUnit = template@spaceUnit, timeUnit = template@timeUnit)

## --- Binary operations -----------------------------------------------------

sparseBinary <- function (op, e1, e2)
{
    first <- isSparseImage(e1)
    second <- isSparseImage(e2)

    if (first && second)
        return(sparseWithSparse(op, e1, e2))
    if (first)
        return(sparseWithOther(op, e1, e2, sparseFirst = TRUE))
    sparseWithOther(op, e2, e1, sparseFirst = FALSE)
}

sparseWithOther <- function (op, x, other, sparseFirst)
{
    ## Anything but a single value has to be matched up position by position,
    ## so there is nothing to be gained by staying packed
    if (!is.atomic(other) || length(other) != 1L)
    {
        dense <- as.array(x)
        return(denseFrom(x, if (sparseFirst) op(dense, asComparable(other)) else op(asComparable(other), dense)))
    }

    zero <- zeroOf(x@values)
    absent <- if (sparseFirst) absentResult(op, zero, other) else absentResult(op, other, zero)

    if (!isStillAbsent(absent))
    {
        dense <- as.array(x)
        return(denseFrom(x, if (sparseFirst) op(dense, other) else op(other, dense)))
    }

    values <- if (sparseFirst) op(x@values, other) else op(other, x@values)
    rebuild(x, x@mask, values)
}

sparseWithSparse <- function (op, e1, e2)
{
    if (!identical(e1@dims, e2@dims))
        stop("Sparse images must have the same dimensions to be combined")
    if (!identical(e1@spatial, e2@spatial))
        stop("Sparse images must have the same number of spatial dimensions")

    absent <- absentResult(op, zeroOf(e1@values), zeroOf(e2@values))
    if (!isStillAbsent(absent))
        return(denseFrom(e1, op(as.array(e1), as.array(e2))))

    ## A result can only be non-zero where at least one operand holds
    ## something, so the union bounds it. Where the operation turns out to
    ## annihilate -- multiplication, say -- tightening afterwards recovers the
    ## intersection without it having to be special-cased here
    combined <- maskCombine(e1@mask, e2@mask, "union")
    locations <- locationCount(e1)
    elements <- elementCount(e1)

    first <- repackValues(e1@mask, e1@values, combined, locations, elements)
    second <- repackValues(e2@mask, e2@values, combined, locations, elements)

    rebuild(e1, combined, op(first, second))
}

## Comparisons against a sparse image produce logical images, which are still
## images; this keeps the storage mode of the operand rather than the operator
asComparable <- function (x) if (isDenseImage(x)) as.array(x) else x

## --- Unary and summary operations ------------------------------------------

sparseMath <- function (fn, x, ...)
{
    zero <- zeroOf(x@values)
    absent <- tryCatch(fn(zero, ...), error = function (e) NULL)

    if (!isStillAbsent(absent))
        return(denseFrom(x, fn(as.array(x), ...)))

    rebuild(x, x@mask, fn(x@values, ...))
}

## Summaries can be answered from the stored values plus a single zero
## standing for every absent location, which is what makes them cheap here.
## tractor.base does the same for its Summary methods, and it is the one place
## its sparse handling is not undone by densification
sparseSummary <- function (fn, x, ..., na.rm = FALSE)
{
    absent <- locationCount(x) > maskCount(x@mask, locationCount(x))
    values <- x@values

    if (absent)
        values <- c(values, zeroOf(values))

    fn(values, ..., na.rm = na.rm)
}

## --- Registration ----------------------------------------------------------

## Registering by hand would mean fifteen operators across three signatures,
## plus the unary and summary functions. method<- is an ordinary function, so
## the registrations are generated instead, each handler closing over the name
## of the function it stands for
registerSparseMethods <- function ()
{
    for (name in binaryOperators)
    {
        generic <- get(name, baseenv())
        ## Only the binary forms are reachable. R resolves a unary minus
        ## through S3 group dispatch, which finds S7's own Ops.S7_object
        ## before any method registered here, so `-x` is not available; write
        ## `0 - x` or `x * -1`, both of which stay sparse
        handler <- local({
            op <- generic
            function (e1, e2) sparseBinary(op, e1, e2)
        })

        S7::`method<-`(generic, list(sparseImage, sparseImage), handler)
        S7::`method<-`(generic, list(sparseImage, S7::class_numeric), handler)
        S7::`method<-`(generic, list(S7::class_numeric, sparseImage), handler)
        S7::`method<-`(generic, list(sparseImage, S7::class_logical), handler)
        S7::`method<-`(generic, list(S7::class_logical, sparseImage), handler)
        S7::`method<-`(generic, list(sparseImage, denseImage), handler)
        S7::`method<-`(generic, list(denseImage, sparseImage), handler)
    }

    for (name in mathFunctions)
    {
        generic <- get(name, baseenv())
        handler <- local({
            fn <- generic
            function (x, ...) sparseMath(fn, x, ...)
        })
        S7::`method<-`(generic, sparseImage, handler)
    }

    for (name in summaryFunctions)
    {
        generic <- get(name, baseenv())
        handler <- local({
            fn <- generic
            function (x, ..., na.rm = FALSE) sparseSummary(fn, x, ..., na.rm = na.rm)
        })
        S7::`method<-`(generic, sparseImage, handler)
    }

    ## NB: `!` cannot be registered, because S7 declines to attach a method to
    ## a primitive that is not an S3 generic in its own right. Use `x == 0`,
    ## or negate the materialised array

    S7::`method<-`(base::mean, sparseImage,
               function (x, ...) sum(x, ...) / length(x))
    S7::`method<-`(base::is.na, sparseImage,
               function (x) array(is.na(as.array(x)), x@dims))

    invisible(NULL)
}
