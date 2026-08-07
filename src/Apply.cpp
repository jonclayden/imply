#include <Rcpp.h>

#include <memory>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Blocks.h"
#include "imply/Sink.h"
#include "imply/Sparse.h"
#include "imply/Narrow.h"

using namespace imply;

namespace {

// How often to look for a pending interrupt. Frequent enough that a slow run
// stops promptly, rare enough that the check itself never shows up in a
// profile: each one is a setjmp and a call, so a hundred calls apart costs
// microseconds over a run of any length
const R_xlen_t interruptInterval = 100;

// Apply an R function over the margins of an array.
//
// Unlike base::apply(), which permutes the whole array into a fresh copy
// before looping, each sub-array is gathered directly through the stride
// vector. Peak memory is therefore the input plus the result, rather than
// twice the input plus the result.
template <typename Accessor, typename Tag>
Rcpp::List applyImpl (const Accessor &source,
                      const std::vector<Extent> &dims, const std::vector<int> &margin,
                      SEXP fun, SEXP callNames, const bool simplify,
                      const R_xlen_t from, const R_xlen_t to,
                      SEXP progress, const R_xlen_t reportEvery, Tag)
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

    // Margin dimensions are taken in the order given, since that is the order
    // the result's dimensions will be in
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

    OffsetWalker marginWalker(marginDims, marginStrides);
    OffsetWalker callWalker(callDims, callStrides);

    // A worker may be given only part of the call space. Seeking straight to
    // its first call avoids walking everything before it
    const R_xlen_t begin = std::max<R_xlen_t>(0, from);
    const R_xlen_t end = (to < 0 ? static_cast<R_xlen_t>(marginWalker.size())
                                 : std::min<R_xlen_t>(to, static_cast<R_xlen_t>(marginWalker.size())));
    const R_xlen_t nCalls = std::max<R_xlen_t>(0, end - begin);
    marginWalker.seek(static_cast<Extent>(begin));

    const R_xlen_t subSize = static_cast<R_xlen_t>(callWalker.size());

    // A sub-array is passed as a bare vector when it has fewer than two
    // dimensions, matching base::apply(). In that case any dimension names
    // become the vector's names instead
    const bool subHasDim = (callDims.size() >= 2);
    Rcpp::RObject subDim, subNames;
    if (subHasDim)
    {
        Rcpp::IntegerVector value(callDims.size());
        for (std::size_t i=0; i<callDims.size(); i++)
            value[i] = static_cast<int>(callDims[i]);
        subDim = value;
        subNames = callNames;
    }
    else if (!Rf_isNull(callNames) && Rf_xlength(callNames) == 1)
        subNames = VECTOR_ELT(callNames, 0);

    // The call is built once and its argument replaced each time round, rather
    // than a fresh call being constructed per iteration
    Rcpp::RObject call = Rf_lang2(fun, R_NilValue);

    // Progress is reported by calling back into R. That is safe here because
    // this loop always runs on the main thread: under forked parallelism each
    // worker is given no reporter, and the parent reports between batches
    // instead. Positions are absolute, so a worker's range still makes sense
    const bool reporting = (!Rf_isNull(progress) && reportEvery > 0);
    Rcpp::RObject progressCall;
    if (reporting)
        progressCall = Rf_lang2(progress, R_NilValue);

    std::unique_ptr<VectorSink> fast;
    std::unique_ptr<ListSink> general;

    for (R_xlen_t k=0; k<nCalls; k++)
    {
        Rcpp::Vector<Tag::sexpType> sub(subSize);
        gather(source, marginWalker.offset(), callWalker, sub.begin());
        if (subHasDim)
        {
            sub.attr("dim") = subDim;
            if (!subNames.isNULL())
                sub.attr("dimnames") = subNames;
        }
        else if (!subNames.isNULL())
            sub.attr("names") = subNames;

        SETCADR(call, sub);
        Rcpp::RObject value = Rf_eval(call, R_GlobalEnv);

        if (general == nullptr && fast == nullptr)
        {
            // The first result decides whether the fast path is available
            const SEXPTYPE type = TYPEOF(value);
            const bool usable = simplify
                                && (type == LGLSXP || type == INTSXP || type == REALSXP || type == CPLXSXP)
                                && Rf_xlength(value) > 0 && !hasAttributes(value);
            if (usable)
                fast = makeVectorSink(type, Rf_xlength(value), nCalls);
            if (fast == nullptr)
                general.reset(new ListSink(nCalls));
        }

        if (fast != nullptr)
        {
            if (!fast->write(k, value))
            {
                // Something no longer fits the preallocated vector, so move
                // what has been written into a list and carry on there
                general.reset(new ListSink(nCalls));
                for (R_xlen_t j=0; j<k; j++)
                    general->write(j, fast->element(j));
                general->write(k, value);
                fast.reset();
            }
        }
        else
            general->write(k, value);

        marginWalker.next();

        if (reporting && ((k + 1) % reportEvery == 0))
        {
            SETCADR(progressCall, Rf_ScalarReal(static_cast<double>(begin + k + 1)));
            Rf_eval(progressCall, R_GlobalEnv);
        }

        // Rcpp runs R_CheckUserInterrupt() inside R_ToplevelExec, so R's
        // longjmp is contained there and an ordinary C++ exception is thrown
        // instead. The stack unwinds properly and the sinks, walkers and
        // buffers below are all destroyed
        if ((k + 1) % interruptInterval == 0)
            Rcpp::checkUserInterrupt();
    }

    if (reporting)
    {
        SETCADR(progressCall, Rf_ScalarReal(static_cast<double>(begin + nCalls)));
        Rf_eval(progressCall, R_GlobalEnv);
    }

    if (fast != nullptr)
        return Rcpp::List::create(Rcpp::Named("values") = fast->finish(),
                                  Rcpp::Named("elementLength") = static_cast<double>(fast->length()),
                                  Rcpp::Named("isList") = false);

    if (general == nullptr)
        general.reset(new ListSink(nCalls));

    return Rcpp::List::create(Rcpp::Named("values") = general->finish(),
                              Rcpp::Named("elementLength") = NA_REAL,
                              Rcpp::Named("isList") = true);
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

std::vector<Extent> dimsFrom (Rcpp::IntegerVector dim)
{
    return std::vector<Extent>(dim.begin(), dim.end());
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List applyOverMargin (Rcpp::RObject x, Rcpp::IntegerVector margin, Rcpp::Function fun,
                            Rcpp::Nullable<Rcpp::List> callNames = R_NilValue, bool simplify = true,
                            double from = 0, double to = -1,
                            Rcpp::Nullable<Rcpp::Function> progress = R_NilValue,
                            double reportEvery = 0)
{
    const std::vector<Extent> dims = dimsOf(x);
    checkLength(x, dims);
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    SEXP names = (callNames.isNull() ? R_NilValue : SEXP(callNames.get()));
    SEXP reporter = (progress.isNull() ? R_NilValue : SEXP(progress.get()));

    Rcpp::List result;
    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;
        result = applyImpl(DenseAccessor<typename Tag::Type>(data), dims, margin0, fun, names, simplify,
                           static_cast<R_xlen_t>(from), static_cast<R_xlen_t>(to),
                           reporter, static_cast<R_xlen_t>(reportEvery), tag);
        return R_NilValue;
    });

    return result;
}

// The same loop over a packed image. Values are widened to double during the
// gather, so the function sees ordinary numbers and never learns that the
// image was stored narrowly
// [[Rcpp::export]]
Rcpp::List applyOverMarginPacked (Rcpp::RawVector values, std::string type, Rcpp::IntegerVector dim,
                                  Rcpp::IntegerVector margin, Rcpp::Function fun,
                                  double slope = 1, double intercept = 0,
                                  Rcpp::Nullable<Rcpp::List> callNames = R_NilValue,
                                  bool simplify = true, double from = 0, double to = -1,
                                  Rcpp::Nullable<Rcpp::Function> progress = R_NilValue,
                                  double reportEvery = 0)
{
    const std::vector<Extent> dims = dimsFrom(dim);
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    SEXP names = (callNames.isNull() ? R_NilValue : SEXP(callNames.get()));
    SEXP reporter = (progress.isNull() ? R_NilValue : SEXP(progress.get()));

    Rcpp::List result;
    dispatchNarrowType(narrowTypeFromName(type), [&](auto stored) -> SEXP {
        typedef decltype(stored) Stored;
        result = applyImpl(NarrowAccessor<Stored>(values.begin(), slope, intercept),
                           dims, margin0, fun, names, simplify,
                           static_cast<R_xlen_t>(from), static_cast<R_xlen_t>(to),
                           reporter, static_cast<R_xlen_t>(reportEvery), RealTag());
        return R_NilValue;
    });

    return result;
}

// ...and over a sparse image, where the gather turns an absent location into
// a zero. Nothing is materialised beyond one sub-array at a time
// [[Rcpp::export]]
Rcpp::List applyOverMarginSparse (Rcpp::RawVector mask, Rcpp::RObject values, Rcpp::IntegerVector dim,
                                  int spatial, Rcpp::IntegerVector margin, Rcpp::Function fun,
                                  Rcpp::Nullable<Rcpp::List> callNames = R_NilValue,
                                  bool simplify = true, double from = 0, double to = -1,
                                  Rcpp::Nullable<Rcpp::Function> progress = R_NilValue,
                                  double reportEvery = 0)
{
    const std::vector<Extent> dims = dimsFrom(dim);
    const std::vector<int> margin0 = checkMargin(margin, static_cast<int>(dims.size()));
    SEXP names = (callNames.isNull() ? R_NilValue : SEXP(callNames.get()));
    SEXP reporter = (progress.isNull() ? R_NilValue : SEXP(progress.get()));

    Extent locations = 1, elements = 1;
    for (int i=0; i<spatial; i++)
        locations *= dims[i];
    for (std::size_t i=spatial; i<dims.size(); i++)
        elements *= dims[i];

    const LocationMask bits(mask, locations);

    Rcpp::List result;
    dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;
        result = applyImpl(SparseAccessor<typename Tag::Type>(bits, packed, elements),
                           dims, margin0, fun, names, simplify,
                           static_cast<R_xlen_t>(from), static_cast<R_xlen_t>(to),
                           reporter, static_cast<R_xlen_t>(reportEvery), tag);
        return R_NilValue;
    });

    return result;
}
