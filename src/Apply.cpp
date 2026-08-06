#include <Rcpp.h>

#include <memory>

#include "Raster.h"
#include "Storage.h"
#include "Dispatch.h"
#include "RImage.h"
#include "Blocks.h"
#include "Sink.h"

using namespace imply;

namespace {

// Apply an R function over the margins of an array.
//
// Unlike base::apply(), which permutes the whole array into a fresh copy
// before looping, each sub-array is gathered directly through the stride
// vector. Peak memory is therefore the input plus the result, rather than
// twice the input plus the result.
template <typename Tag>
Rcpp::List applyImpl (SEXP x, const typename Tag::type *data,
                      const std::vector<Extent> &dims, const std::vector<int> &margin,
                      SEXP fun, SEXP callNames, const bool simplify, Tag)
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

    offsetWalker marginWalker(marginDims, marginStrides);
    offsetWalker callWalker(callDims, callStrides);

    const R_xlen_t nCalls = static_cast<R_xlen_t>(marginWalker.size());
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

    std::unique_ptr<vectorSink> fast;
    std::unique_ptr<listSink> general;

    for (R_xlen_t k=0; k<nCalls; k++)
    {
        Rcpp::Vector<Tag::sexpType> sub(subSize);
        gather(data, marginWalker.offset(), callWalker, sub.begin());
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
                general.reset(new listSink(nCalls));
        }

        if (fast != nullptr)
        {
            if (!fast->write(k, value))
            {
                // Something no longer fits the preallocated vector, so move
                // what has been written into a list and carry on there
                general.reset(new listSink(nCalls));
                for (R_xlen_t j=0; j<k; j++)
                    general->write(j, fast->element(j));
                general->write(k, value);
                fast.reset();
            }
        }
        else
            general->write(k, value);

        marginWalker.next();
    }

    if (fast != nullptr)
        return Rcpp::List::create(Rcpp::Named("values") = fast->finish(),
                                  Rcpp::Named("elementLength") = static_cast<double>(fast->length()),
                                  Rcpp::Named("isList") = false);

    if (general == nullptr)
        general.reset(new listSink(nCalls));

    return Rcpp::List::create(Rcpp::Named("values") = general->finish(),
                              Rcpp::Named("elementLength") = NA_REAL,
                              Rcpp::Named("isList") = true);
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List applyOverMargin (Rcpp::RObject x, Rcpp::IntegerVector margin, Rcpp::Function fun,
                            Rcpp::Nullable<Rcpp::List> callNames = R_NilValue, bool simplify = true)
{
    const std::vector<Extent> dims = dimsOf(x);
    checkLength(x, dims);
    const int nDims = static_cast<int>(dims.size());

    std::vector<int> margin0;
    margin0.reserve(margin.size());
    for (R_xlen_t i=0; i<margin.size(); i++)
    {
        if (margin[i] == NA_INTEGER || margin[i] < 1 || margin[i] > nDims)
            Rcpp::stop("Margin %d is out of range for an array with %d dimensions", i+1, nDims);
        margin0.push_back(margin[i] - 1);
    }

    for (std::size_t i=0; i<margin0.size(); i++)
    {
        for (std::size_t j=i+1; j<margin0.size(); j++)
        {
            if (margin0[i] == margin0[j])
                Rcpp::stop("Margin contains a repeated dimension");
        }
    }

    SEXP names = (callNames.isNull() ? R_NilValue : SEXP(callNames.get()));

    Rcpp::List result;
    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        result = applyImpl(x, data, dims, margin0, fun, names, simplify, tag);
        return R_NilValue;
    });

    return result;
}
