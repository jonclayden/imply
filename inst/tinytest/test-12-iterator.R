## StrideIterator: a random-access iterator over a strided run of values,
## exercised here via std algorithms (accumulate, reverse, sort) that require
## a real iterator, not just an offset walker.

## These are internal probes for now, exactly like the Raster ones in
## test-01-raster.R
strideIteratorRead <- imply:::strideIteratorRead
strideIteratorAccumulate <- imply:::strideIteratorAccumulate
strideIteratorReverse <- imply:::strideIteratorReverse
strideIteratorSort <- imply:::strideIteratorSort

x <- c(5, 3, 8, 1, 9, 2, 7, 4, 6, 0)

## --- Reading -----------------------------------------------------------

## Stride 1 is a plain, contiguous walk
expect_equal(strideIteratorRead(x, 1), x)

## A stride visits every nth element, same as R's own subsetting
for (stride in 2:4)
    expect_equal(strideIteratorRead(x, stride), x[seq(1, length(x), by=stride)],
                 info=paste("stride", stride))

## A stride that doesn't evenly divide the length still stops at the last
## element it can fully reach, not past it
expect_equal(strideIteratorRead(1:7, 3), c(1,4,7))

## --- Algorithms that need more than a forward walk ----------------------

## std::accumulate only needs input iteration, but exercises operator+ on a
## non-trivial stride
for (stride in 1:3)
    expect_equal(strideIteratorAccumulate(x, stride), sum(x[seq(1, length(x), by=stride)]),
                 info=paste("stride", stride))

## std::reverse needs bidirectional access and mutation through the iterator;
## elements not on the stride are untouched
for (stride in 1:3) {
    idx <- seq(1, length(x), by=stride)
    expected <- x
    expected[idx] <- rev(x[idx])
    expect_equal(strideIteratorReverse(x, stride), expected, info=paste("stride", stride))
}

## std::sort requires a genuine RandomAccessIterator (operator[], operator-,
## the full set of comparisons) -- this is the strongest check that the
## iterator category is honestly random-access, not just claimed to be
for (stride in 1:3) {
    idx <- seq(1, length(x), by=stride)
    expected <- x
    expected[idx] <- sort(x[idx])
    expect_equal(strideIteratorSort(x, stride), expected, info=paste("stride", stride))
}

## A single-element run is a valid (degenerate) case for every algorithm above
expect_equal(strideIteratorRead(x, length(x)), x[1])
expect_equal(strideIteratorSort(x, length(x)), x)
