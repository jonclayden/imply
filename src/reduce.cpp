#include <Rcpp.h>

#include <cmath>
#include <string>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Blocks.h"
#include "imply/Sparse.h"
#include "imply/Narrow.h"
#include "imply/Parallel.h"

using namespace imply;

namespace {

// Reductions that can be computed in one pass, without an R callback. This is
// the tier that genuinely parallelises: an R function cannot be called from a
// worker thread, but none of these touch R at all
enum class Reduction
{
    sum, mean, minimum, maximum, range, product, variance, deviation,
    whichMinimum, whichMaximum, any, all, countNA
};

Reduction reductionFromName (const std::string &name)
{
    if (name == "sum")      return Reduction::sum;
    if (name == "mean")     return Reduction::mean;
    if (name == "min")      return Reduction::minimum;
    if (name == "max")      return Reduction::maximum;
    if (name == "range")    return Reduction::range;
    if (name == "prod")     return Reduction::product;
    if (name == "var")      return Reduction::variance;
    if (name == "sd")       return Reduction::deviation;
    if (name == "which.min") return Reduction::whichMinimum;
    if (name == "which.max") return Reduction::whichMaximum;
    if (name == "any")      return Reduction::any;
    if (name == "all")      return Reduction::all;
    if (name == "countNA")  return Reduction::countNA;
    Rcpp::stop("Unknown reduction \"%s\"", name);
}

// How many values each call contributes to the result
int reductionWidth (const Reduction what)
{
    return (what == Reduction::range ? 2 : 1);
}

// Everything accumulates in double, whatever the values were stored as, so
// the answer does not depend on the storage type and error does not compound
struct Accumulator
{
    double total = 0.0;
    double sumSquares = 0.0;
    double low = R_PosInf;
    double high = R_NegInf;
    double logProduct = 0.0;
    double product = 1.0;
    Extent lowAt = 0, highAt = 0, used = 0, missing = 0;
    bool anyTrue = false, allTrue = true;

    void add (const double value, const Extent index)
    {
        if (ISNAN(value))
        {
            missing++;
            return;
        }

        total += value;
        sumSquares += value * value;
        product *= value;

        if (value < low)  { low = value;  lowAt = index; }
        if (value > high) { high = value; highAt = index; }

        if (value != 0.0) anyTrue = true;
        else              allTrue = false;

        used++;
    }
};

void writeResult (const Accumulator &a, const Reduction what, const bool naRm,
                  double * const out, const R_xlen_t position, const int width)
{
    const bool spoiled = (a.missing > 0 && !naRm);
    const bool empty = (a.used == 0);

    double first = NA_REAL, second = NA_REAL;

    if (what == Reduction::countNA)
        first = static_cast<double>(a.missing);
    else if (spoiled)
        ; // leave as NA
    else
    {
        switch (what)
        {
            case Reduction::sum:      first = a.total; break;
            case Reduction::product:  first = (empty ? 1.0 : a.product); break;
            case Reduction::mean:     first = (empty ? R_NaN : a.total / double(a.used)); break;
            case Reduction::minimum:  first = (empty ? R_PosInf : a.low); break;
            case Reduction::maximum:  first = (empty ? R_NegInf : a.high); break;

            case Reduction::range:
            first = (empty ? R_PosInf : a.low);
            second = (empty ? R_NegInf : a.high);
            break;

            case Reduction::variance:
            case Reduction::deviation:
            {
                if (a.used < 2)
                    first = NA_REAL;
                else
                {
                    const double n = double(a.used);
                    const double mean = a.total / n;
                    // The corrected two-pass form is not available in one
                    // pass, so the sum of squares is used with the mean
                    // subtracted afterwards, clamped at zero against
                    // cancellation
                    double value = (a.sumSquares - n * mean * mean) / (n - 1.0);
                    if (value < 0.0)
                        value = 0.0;
                    first = (what == Reduction::deviation ? std::sqrt(value) : value);
                }
                break;
            }

            // Positions are one-based, as R reports them
            case Reduction::whichMinimum: first = (empty ? NA_REAL : double(a.lowAt) + 1.0); break;
            case Reduction::whichMaximum: first = (empty ? NA_REAL : double(a.highAt) + 1.0); break;

            case Reduction::any: first = (a.anyTrue ? 1.0 : 0.0); break;
            case Reduction::all: first = (a.allTrue ? 1.0 : 0.0); break;
            case Reduction::countNA: break;
        }
    }

    // any() and all() have the three-valued behaviour R gives them: a missing
    // value only matters when it could change the answer
    if (!naRm && a.missing > 0)
    {
        if (what == Reduction::any)
            first = (a.anyTrue ? 1.0 : NA_REAL);
        else if (what == Reduction::all)
            first = (a.allTrue ? NA_REAL : 0.0);
    }

    out[position * width] = first;
    if (width > 1)
        out[position * width + 1] = second;
}

// The same margin decomposition as the apply engine, but with the R callback
// replaced by an accumulator, which is what allows the loop to run on worker
// threads
template <typename Accessor>
void reduceImpl (const Accessor &source, const std::vector<Extent> &dims,
                 const std::vector<int> &margin, const Reduction what, const bool naRm,
                 const int threads, double * const out,
                 const Extent from, const Extent to)
{
    const int nDims = static_cast<int>(dims.size());

    std::vector<Extent> strides(nDims);
    Extent stride = 1;
    for (int i=0; i<nDims; i++)
    {
        strides[i] = stride;
        stride *= dims[i];
    }

    std::vector<bool> retained(nDims, false);
    for (std::size_t i=0; i<margin.size(); i++)
        retained[margin[i]] = true;

    std::vector<Extent> marginDims, marginStrides, callDims, callStrides;
    for (std::size_t i=0; i<margin.size(); i++)
    {
        marginDims.push_back(dims[margin[i]]);
        marginStrides.push_back(strides[margin[i]]);
    }
    for (int i=0; i<nDims; i++)
    {
        if (!retained[i])
        {
            callDims.push_back(dims[i]);
            callStrides.push_back(strides[i]);
        }
    }

    const OffsetWalker marginTemplate(marginDims, marginStrides);
    const OffsetWalker callTemplate(callDims, callStrides);
    const Extent first = from;
    const Extent last = std::min(to, marginTemplate.size());
    const int width = reductionWidth(what);

    if (last <= first)
        return;

    // Each chunk owns a disjoint run of calls and writes only its own slice of
    // the output, so there is nothing to synchronise
    parallelFor(last - first, threads, [&](const Extent begin, const Extent end) {
        OffsetWalker margins = marginTemplate;
        OffsetWalker values = callTemplate;
        margins.seek(first + begin);

        for (Extent k=first+begin; k<first+end; k++)
        {
            const Offset base = margins.offset();
            const Extent n = values.size();

            Accumulator a;
            values.reset();
            for (Extent i=0; i<n; i++)
            {
                a.add(source[static_cast<Extent>(base + values.offset())], i);
                values.next();
            }

            writeResult(a, what, naRm, out, static_cast<R_xlen_t>(k), width);
            margins.next();
        }
    });
}

std::vector<int> checkMargin (Rcpp::IntegerVector margin, const int nDims)
{
    std::vector<int> result;
    result.reserve(margin.size());

    for (R_xlen_t i=0; i<margin.size(); i++)
    {
        if (margin[i] == NA_INTEGER || margin[i] < 1 || margin[i] > nDims)
            Rcpp::stop("Margin %d is out of range for an array with %d dimensions", i+1, nDims);
        result.push_back(margin[i] - 1);
    }

    for (std::size_t i=0; i<result.size(); i++)
    {
        for (std::size_t j=i+1; j<result.size(); j++)
        {
            if (result[i] == result[j])
                Rcpp::stop("Margin contains a repeated dimension");
        }
    }

    return result;
}

Extent callCount (const std::vector<Extent> &dims, const std::vector<int> &margin)
{
    Extent n = 1;
    for (std::size_t i=0; i<margin.size(); i++)
        n *= dims[margin[i]];
    return n;
}

// Nothing inside the kernel may touch the R API, since it runs on worker
// threads, so the interrupt is looked for between slabs of calls instead.
// Slabbing costs one parallel launch per slab, which is only worth paying when
// there is enough work for it to disappear into the noise
template <typename Accessor>
void runReduction (const Accessor &source, const std::vector<Extent> &dims,
                   const std::vector<int> &margin, const Reduction what, const bool naRm,
                   const int threads, double * const out)
{
    const Extent nCalls = callCount(dims, margin);
    Extent total = 1;
    for (std::size_t i=0; i<dims.size(); i++)
        total *= dims[i];

    const Extent slabs = (total > 10000000u ? std::min<Extent>(nCalls, 32) : 1);

    if (slabs <= 1)
    {
        reduceImpl(source, dims, margin, what, naRm, threads, out, 0, nCalls);
        return;
    }

    const Extent size = (nCalls + slabs - 1) / slabs;
    for (Extent from=0; from<nCalls; from+=size)
    {
        reduceImpl(source, dims, margin, what, naRm, threads, out, from, std::min(from + size, nCalls));
        Rcpp::checkUserInterrupt();
    }
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::NumericVector reduceOverMargin (Rcpp::RObject x, Rcpp::IntegerVector margin, std::string what,
                                      bool naRm = false, int threads = 0)
{
    const std::vector<Extent> dims = dimsOf(x);
    checkLength(x, dims);
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    const Reduction kind = reductionFromName(what);

    Rcpp::NumericVector result(static_cast<R_xlen_t>(callCount(dims, margin0)) * reductionWidth(kind));

    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;
        if constexpr (Tag::kind == StorageType::complex)
            Rcpp::stop("Complex data are not supported by imreduce()");
        else
            runReduction(DenseAccessor<typename Tag::Type>(data), dims, margin0, kind, naRm,
                         threads, result.begin());
        return R_NilValue;
    });

    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduceOverMarginPacked (Rcpp::RawVector values, std::string type,
                                            Rcpp::IntegerVector dim, Rcpp::IntegerVector margin,
                                            std::string what, double slope = 1, double intercept = 0,
                                            bool naRm = false, int threads = 0)
{
    const std::vector<Extent> dims(dim.begin(), dim.end());
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    const Reduction kind = reductionFromName(what);

    Rcpp::NumericVector result(static_cast<R_xlen_t>(callCount(dims, margin0)) * reductionWidth(kind));

    dispatchNarrowType(narrowTypeFromName(type), [&](auto stored) -> SEXP {
        typedef decltype(stored) Stored;
        runReduction(NarrowAccessor<Stored>(values.begin(), slope, intercept), dims, margin0, kind,
                     naRm, threads, result.begin());
        return R_NilValue;
    });

    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector reduceOverMarginSparse (Rcpp::RawVector mask, Rcpp::RObject values,
                                            Rcpp::IntegerVector dim, int spatial,
                                            Rcpp::IntegerVector margin, std::string what,
                                            bool naRm = false, int threads = 0)
{
    const std::vector<Extent> dims(dim.begin(), dim.end());
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    const Reduction kind = reductionFromName(what);

    Extent locations = 1, elements = 1;
    for (int i=0; i<spatial; i++)
        locations *= dims[i];
    for (std::size_t i=spatial; i<dims.size(); i++)
        elements *= dims[i];

    const LocationMask bits(mask, locations);
    Rcpp::NumericVector result(static_cast<R_xlen_t>(callCount(dims, margin0)) * reductionWidth(kind));

    dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;
        if constexpr (Tag::kind == StorageType::complex)
            Rcpp::stop("Complex data are not supported by imreduce()");
        else
            runReduction(SparseAccessor<typename Tag::Type>(bits, packed, elements), dims, margin0,
                         kind, naRm, threads, result.begin());
        return R_NilValue;
    });

    return result;
}
