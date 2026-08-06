#include <Rcpp.h>

// Supplies Rcpp::wrap() for std::array, so a fixed-dimensionality raster's
// extents convert to R exactly as a dynamic one's std::vector does
#include "RcppArray.h"

#include "Raster.h"
#include "Storage.h"
#include "Dispatch.h"
#include "RImage.h"
#include "Parallel.h"

using namespace imply;

namespace {

template <typename Raster>
Rcpp::List rasterInfoImpl (const Raster &r)
{
    return Rcpp::List::create(
        Rcpp::Named("dim") = Rcpp::wrap(r.dim()),
        Rcpp::Named("strides") = Rcpp::wrap(r.strides()),
        Rcpp::Named("nDims") = r.nDims(),
        Rcpp::Named("spatial") = r.spatial(),
        Rcpp::Named("size") = static_cast<double>(r.size()),
        Rcpp::Named("spatialSize") = static_cast<double>(r.spatialSize()),
        Rcpp::Named("elementSize") = static_cast<double>(r.elementSize()),
        Rcpp::Named("contiguous") = r.isContiguous(),
        Rcpp::Named("fixed") = Raster::isFixed);
}

// Sum along every line running in the given direction. Lines are enumerated
// with the first remaining dimension moving fastest, which matches the order
// base::apply() produces over the complementary margins
template <typename Raster, typename Tag>
SEXP lineSumsImpl (const Raster &r, const typename Tag::type *data, const int dim, Tag)
{
    if constexpr (Tag::kind == storageType::complex)
    {
        Rcpp::stop("Complex data are not yet supported by lineSums()");
        return R_NilValue;
    }
    else
    {
        const Extent nLines = r.countLines(dim);
        const Extent n = r.dim(dim);
        const Extent stride = r.stride(dim);

        Rcpp::NumericVector result(nLines);
        for (Extent i=0; i<nLines; i++)
        {
            const Offset base = r.lineOffset(i, dim);
            double sum = 0.0;
            bool missing = false;

            for (Extent j=0; j<n; j++)
            {
                const typename Tag::type value = data[base + static_cast<Offset>(j*stride)];
                if (Tag::isNA(value))
                {
                    missing = true;
                    break;
                }
                sum += static_cast<double>(value);
            }

            result[i] = missing ? NA_REAL : sum;
        }

        return result;
    }
}

// Materialise a permuted view. Nothing is permuted in memory: the view carries
// reordered strides, and the walk below reads through them
template <typename Raster, typename Tag>
SEXP permuteImpl (const Raster &r, const typename Tag::type *data, const std::vector<int> &order,
                  const int threads, Tag)
{
    const Raster permuted = r.permute(order);

    // Allocated here, on the main thread, because nothing inside the parallel
    // region below may touch the R API
    Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(permuted.size()));
    typename Tag::type * const out = result.begin();

    // Each chunk owns a disjoint range of the output, so there is nothing to
    // synchronise and the result does not depend on how the work is divided
    parallelFor(permuted.size(), threads, [&](const Extent begin, const Extent end) {
        // Declared inside, so each worker has its own. Hoisting it out of the
        // inner loop still matters: building the index per element would cost
        // an allocation per element on the runtime-dimensionality path
        typename Raster::index loc;
        if constexpr (!Raster::isFixed)
            loc.resize(permuted.nDims());

        for (Extent n=begin; n<end; n++)
        {
            permuted.expandIndex(n, loc);
            out[n] = data[permuted.flattenIndex(loc)];
        }
    });

    result.attr("dim") = Rcpp::wrap(permuted.dim());
    return result;
}

} // anonymous namespace

// [[Rcpp::export]]
SEXP rasterInfo (Rcpp::RObject x, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, bool forceDynamic = false)
{
    const rasterSpec spec = specOf(x, spatial);

    if (forceDynamic)
        return rasterInfoImpl(dynamicRaster(spec.dims, spec.spatial));

    return dispatchDims(spec.nDims(), [&](auto tag) -> SEXP {
        return rasterInfoImpl(raster<decltype(tag)::value>(spec.dims, spec.spatial));
    });
}

// Convert one-based array locations (a matrix, one row per location) to
// one-based linear indices
// [[Rcpp::export]]
SEXP flattenIndices (Rcpp::RObject x, Rcpp::IntegerMatrix locs, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, bool forceDynamic = false)
{
    const rasterSpec spec = specOf(x, spatial);
    const int nDims = spec.nDims();

    if (locs.ncol() != nDims)
        Rcpp::stop("Location matrix has %d columns, but the object has %d dimensions", locs.ncol(), nDims);

    const R_xlen_t nLocs = locs.nrow();
    Rcpp::NumericVector result(nLocs);

    auto run = [&](const auto &r) {
        typename std::decay_t<decltype(r)>::index loc;
        if constexpr (!std::decay_t<decltype(r)>::isFixed)
            loc.resize(nDims);

        for (R_xlen_t i=0; i<nLocs; i++)
        {
            for (int j=0; j<nDims; j++)
            {
                const int value = locs(i,j);
                if (value == NA_INTEGER || value < 1 || static_cast<Extent>(value) > spec.dims[j])
                    Rcpp::stop("Location [%d,%d] is out of range", i+1, j+1);
                loc[j] = static_cast<Extent>(value - 1);
            }
            result[i] = static_cast<double>(r.flattenIndex(loc)) + 1.0;
        }
    };

    if (forceDynamic)
        run(dynamicRaster(spec.dims, spec.spatial));
    else
        dispatchDims(nDims, [&](auto tag) -> SEXP {
            run(raster<decltype(tag)::value>(spec.dims, spec.spatial));
            return R_NilValue;
        });

    return result;
}

// The inverse: one-based linear indices to a matrix of one-based locations
// [[Rcpp::export]]
SEXP expandIndices (Rcpp::RObject x, Rcpp::NumericVector indices, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, bool forceDynamic = false)
{
    const rasterSpec spec = specOf(x, spatial);
    const int nDims = spec.nDims();

    const R_xlen_t n = indices.size();
    Rcpp::IntegerMatrix result(n, nDims);

    auto run = [&](const auto &r) {
        typename std::decay_t<decltype(r)>::index loc;
        if constexpr (!std::decay_t<decltype(r)>::isFixed)
            loc.resize(nDims);

        for (R_xlen_t i=0; i<n; i++)
        {
            if (Rcpp::NumericVector::is_na(indices[i]) || indices[i] < 1 || indices[i] > static_cast<double>(r.size()))
                Rcpp::stop("Index %d is out of range", i+1);
            r.expandIndex(static_cast<Extent>(indices[i]) - 1, loc);
            for (int j=0; j<nDims; j++)
                result(i,j) = static_cast<int>(loc[j]) + 1;
        }
    };

    if (forceDynamic)
        run(dynamicRaster(spec.dims, spec.spatial));
    else
        dispatchDims(nDims, [&](auto tag) -> SEXP {
            run(raster<decltype(tag)::value>(spec.dims, spec.spatial));
            return R_NilValue;
        });

    return result;
}

// [[Rcpp::export]]
SEXP lineSums (Rcpp::RObject x, int dim, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, bool forceDynamic = false)
{
    const rasterSpec spec = specOf(x, spatial);

    if (dim < 1 || dim > spec.nDims())
        Rcpp::stop("Dimension %d is out of range", dim);
    const int dim0 = dim - 1;

    return dispatchType(x, [&](auto typeTag, auto *data) -> SEXP {
        if (forceDynamic)
            return lineSumsImpl(dynamicRaster(spec.dims, spec.spatial), data, dim0, typeTag);

        return dispatchDims(spec.nDims(), [&](auto dimTag) -> SEXP {
            return lineSumsImpl(raster<decltype(dimTag)::value>(spec.dims, spec.spatial), data, dim0, typeTag);
        });
    });
}

// Partition the spatial locations into runs, either sized to a target number of
// values or split into a fixed number of pieces
// [[Rcpp::export]]
SEXP blockPartition (Rcpp::RObject x, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, double targetElements = 65536, int count = 0)
{
    const rasterSpec spec = specOf(x, spatial);
    const dynamicRaster r(spec.dims, spec.spatial);

    const std::vector<block> blocks = (count > 0
        ? r.blocksForCount(static_cast<Extent>(count))
        : r.blocks(static_cast<Extent>(targetElements)));

    Rcpp::NumericMatrix result(static_cast<int>(blocks.size()), 2);
    Rcpp::colnames(result) = Rcpp::CharacterVector::create("start", "length");
    for (std::size_t i=0; i<blocks.size(); i++)
    {
        result(static_cast<int>(i), 0) = static_cast<double>(blocks[i].start) + 1.0;
        result(static_cast<int>(i), 1) = static_cast<double>(blocks[i].length);
    }

    return result;
}

// Which parallel backend was compiled in, for the test suite and for
// reporting to the user
// [[Rcpp::export]]
Rcpp::List parallelInfo ()
{
    return Rcpp::List::create(
        Rcpp::Named("backend") = std::string(parallelBackend()),
        Rcpp::Named("available") = parallelAvailable());
}

// [[Rcpp::export]]
SEXP chunkPartition (double items, int threads)
{
    return Rcpp::wrap(static_cast<double>(chunkCount(static_cast<Extent>(items), threads)));
}

// [[Rcpp::export]]
SEXP permuteView (Rcpp::RObject x, Rcpp::IntegerVector order, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue, bool forceDynamic = false, int threads = 0)
{
    const rasterSpec spec = specOf(x, spatial);
    const int nDims = spec.nDims();

    if (order.size() != nDims)
        Rcpp::stop("Permutation has length %d, but the object has %d dimensions", order.size(), nDims);

    std::vector<int> order0(nDims);
    for (int i=0; i<nDims; i++)
    {
        if (order[i] == NA_INTEGER || order[i] < 1 || order[i] > nDims)
            Rcpp::stop("Permutation entry %d is out of range", i+1);
        order0[i] = order[i] - 1;
    }

    return dispatchType(x, [&](auto typeTag, auto *data) -> SEXP {
        if (forceDynamic)
            return permuteImpl(dynamicRaster(spec.dims, spec.spatial), data, order0, threads, typeTag);

        return dispatchDims(nDims, [&](auto dimTag) -> SEXP {
            return permuteImpl(raster<decltype(dimTag)::value>(spec.dims, spec.spatial), data, order0, threads, typeTag);
        });
    });
}

// Address of the object's data block, used by the test suite to confirm that
// read-only paths borrow R's memory rather than duplicating it
// [[Rcpp::export]]
std::string dataAddress (Rcpp::RObject x)
{
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%p", DATAPTR_RO(x));
    return std::string(buffer);
}
