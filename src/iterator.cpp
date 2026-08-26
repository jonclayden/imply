#include <Rcpp.h>

#include <algorithm>
#include <numeric>

#include "imply/Iterator.h"

using namespace imply;

// Internal probes for StrideIterator, exercised from test-12-iterator.R.
// Each takes a stride in R's 1-based counting, so stride=1 is every element,
// stride=2 is every other, and so on -- the same convention as R's own
// seq(1, length(x), by=stride)

namespace {

R_xlen_t strideLength (const R_xlen_t n, const int stride)
{
    return (n + stride - 1) / stride;
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::NumericVector strideIteratorRead (Rcpp::NumericVector x, int stride)
{
    const R_xlen_t n = strideLength(x.length(), stride);
    StrideIterator<double> it(REAL(x), stride);

    Rcpp::NumericVector result(n);
    for (R_xlen_t i=0; i<n; i++, ++it)
        result[i] = *it;
    return result;
}

// [[Rcpp::export]]
double strideIteratorAccumulate (Rcpp::NumericVector x, int stride)
{
    const R_xlen_t n = strideLength(x.length(), stride);
    StrideIterator<double> begin(REAL(x), stride);
    return std::accumulate(begin, begin + n, 0.0);
}

// [[Rcpp::export]]
Rcpp::NumericVector strideIteratorReverse (Rcpp::NumericVector x, int stride)
{
    Rcpp::NumericVector result = Rcpp::clone(x);
    const R_xlen_t n = strideLength(result.length(), stride);
    StrideIterator<double> begin(REAL(result), stride);
    std::reverse(begin, begin + n);
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector strideIteratorSort (Rcpp::NumericVector x, int stride)
{
    Rcpp::NumericVector result = Rcpp::clone(x);
    const R_xlen_t n = strideLength(result.length(), stride);
    StrideIterator<double> begin(REAL(result), stride);
    std::sort(begin, begin + n);
    return result;
}
