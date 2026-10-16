#ifndef _IMPLY_RIMAGE_H_
#define _IMPLY_RIMAGE_H_

#include <Rcpp.h>

#include <vector>

#include "Raster.h"

namespace imply {

// Read a dimension vector from an R object. A plain vector with no dim
// attribute is treated as one-dimensional, so that no custom class is ever
// required of the caller.
//
// Constructing an IntegerVector from the attribute coerces a double-valued dim
// (which R itself will accept) and protects the intermediate, both of which
// have to be done by hand through the C API
inline std::vector<Extent> dimsOf (const Rcpp::RObject &x)
{
    if (!x.hasAttribute("dim"))
        return std::vector<Extent>(1, static_cast<Extent>(Rf_xlength(x)));

    const Rcpp::IntegerVector dim(x.attr("dim"));
    std::vector<Extent> result;
    result.reserve(dim.size());

    for (R_xlen_t i=0; i<dim.size(); i++)
    {
        if (dim[i] == NA_INTEGER || dim[i] < 0)
            Rcpp::stop("Dimensions must not be missing or negative");
        result.push_back(static_cast<Extent>(dim[i]));
    }

    return result;
}

// Resolve the spatial/element split. An explicit argument wins; otherwise the
// grid of the image's geometry is used; otherwise the leading three dimensions
// (or all of them, if fewer) are taken to be spatial.
//
// S7 stores properties as ordinary attributes, so the geometry is an attribute
// of the image and its grid an attribute of that, and no S7 knowledge is
// needed here. Neither lookup allocates, so nothing needs protecting
inline int spatialOf (const Rcpp::RObject &x, const Rcpp::Nullable<Rcpp::IntegerVector> &spatial, const int nDims)
{
    int result = NA_INTEGER;

    if (spatial.isNotNull())
    {
        const Rcpp::IntegerVector value(spatial.get());
        if (value.size() != 1)
            Rcpp::stop("Number of spatial dimensions must be a single value");
        result = value[0];
    }
    else
    {
        const SEXP geometry = Rf_getAttrib(x, Rf_install("geometry"));
        if (geometry != R_NilValue)
        {
            const SEXP grid = Rf_getAttrib(geometry, Rf_install("dims"));
            if (grid != R_NilValue)
                result = static_cast<int>(Rf_xlength(grid));
        }
    }

    if (result == NA_INTEGER || result < 0)
        result = std::min(3, nDims);

    if (result > nDims)
        Rcpp::stop("Number of spatial dimensions (%d) exceeds the dimensionality (%d)", result, nDims);

    return result;
}

// Total element count implied by a dimension vector, checked against the
// object's actual length
inline void checkLength (const Rcpp::RObject &x, const std::vector<Extent> &dims)
{
    Extent length = 1;
    for (std::size_t i=0; i<dims.size(); i++)
        length *= dims[i];

    const Extent actual = static_cast<Extent>(Rf_xlength(x));
    if (length != actual)
        Rcpp::stop("Dimensions imply %d elements, but the object has %d", length, actual);
}

// The three steps above are always taken together, so bundle them
struct RasterSpec
{
    std::vector<Extent> dims;
    int spatial;

    int nDims () const { return static_cast<int>(dims.size()); }
};

inline RasterSpec specOf (const Rcpp::RObject &x, const Rcpp::Nullable<Rcpp::IntegerVector> &spatial)
{
    RasterSpec spec;
    spec.dims = dimsOf(x);
    checkLength(x, spec.dims);
    spec.spatial = spatialOf(x, spatial, spec.nDims());
    return spec;
}

} // namespace imply

#endif
