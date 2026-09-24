#include <Rcpp.h>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/Sparse.h"
#include "imply/View.h"

using namespace imply;

// Rearrange a vector held in storage order into view order, as an array with
// the view's dimensions. This is how a dense image is reoriented, and how a
// sparse image's mask is shown in the image's own axis order
// [[Rcpp::export]]
SEXP viewGather (Rcpp::RObject x, Rcpp::IntegerVector dim, SEXP layout = R_NilValue)
{
    const ViewMap view(std::vector<Extent>(dim.begin(), dim.end()), layoutFrom(layout));
    if (static_cast<R_xlen_t>(view.size()) != Rf_xlength(x))
        Rcpp::stop("Dimensions imply %d elements, but the object has %d", double(view.size()), double(Rf_xlength(x)));

    return dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;
        Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(view.size()));
        gatherView(DenseAccessor<typename Tag::Type>(data), view, result.begin());
        result.attr("dim") = dim;
        return result;
    });
}

// Map one-based linear indices in view order to one-based linear indices in
// storage order, so that element access need not materialise anything
// [[Rcpp::export]]
Rcpp::NumericVector viewIndices (Rcpp::NumericVector indices, Rcpp::IntegerVector dim, SEXP layout = R_NilValue)
{
    const ViewMap view(std::vector<Extent>(dim.begin(), dim.end()), layoutFrom(layout));
    const double total = static_cast<double>(view.size());
    Rcpp::NumericVector result(indices.size());

    for (R_xlen_t k=0; k<indices.size(); k++)
    {
        const double index = indices[k];
        if (Rcpp::NumericVector::is_na(index) || index < 1 || index > total)
            Rcpp::stop("Index %d is out of range", double(k + 1));
        result[k] = static_cast<double>(view.storageIndex(static_cast<Extent>(index) - 1) + 1);
    }

    return result;
}
