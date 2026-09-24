#ifndef _IMPLY_VIEW_H_
#define _IMPLY_VIEW_H_

#include <Rcpp.h>

#include <cstdlib>
#include <vector>

#include "Raster.h"
#include "Blocks.h"

namespace imply {

// How an image's axes map onto its storage. The image's own axes, the ones
// indexing, dim() and the geometry refer to, are the view; the order the
// values sit in memory (or in a file) is the storage.
//
// A layout gives, for each view axis, the one-based rank of the storage axis
// it corresponds to, negated if the view runs along it backwards. Rank one is
// the fastest-varying storage axis. So c(1, 2, 3) is the identity, c(-1, 2, 3)
// reverses the first axis, and c(2, 3, 4, 1) describes a four-dimensional
// image whose fourth axis is interleaved, varying fastest in storage. An empty
// layout means the identity. The convention is MRtrix's, but one-based, so
// that a reversed first axis can be written down.
//
// Reorienting an image therefore changes only its layout and geometry, and
// every kernel keeps working, because a view is nothing more than signed
// strides and a base offset for the walkers in Blocks.h
class ViewMap
{
protected:
    std::vector<int> ranks_;
    std::vector<bool> reversed_;

public:
    std::vector<Extent> dims;           // view dimensions
    std::vector<Extent> storageDims;    // the same extents, in storage order
    std::vector<Offset> strides;        // storage step for each view axis
    Offset base;                        // storage offset of the view's first element
    std::vector<Offset> inverseStrides; // view step for each storage axis
    Offset inverseBase;                 // view offset of the first stored element
    bool identity;

    ViewMap () : base(0), inverseBase(0), identity(true) {}

    explicit ViewMap (const std::vector<Extent> &dims, const std::vector<int> &layout = std::vector<int>())
        : dims(dims), base(0), inverseBase(0), identity(true)
    {
        const std::size_t n = dims.size();
        ranks_.resize(n);
        reversed_.assign(n, false);

        if (layout.empty())
        {
            for (std::size_t i=0; i<n; i++)
                ranks_[i] = static_cast<int>(i);
        }
        else
        {
            if (layout.size() != n)
                Rcpp::stop("Layout must have one entry per dimension (%d)", static_cast<int>(n));

            std::vector<bool> seen(n, false);
            for (std::size_t i=0; i<n; i++)
            {
                const int value = layout[i];
                if (value == NA_INTEGER || value == 0 || std::abs(value) > static_cast<int>(n))
                    Rcpp::stop("Layout entries must be non-zero storage ranks between 1 and %d", static_cast<int>(n));
                const int rank = std::abs(value) - 1;
                if (seen[rank])
                    Rcpp::stop("Layout must name each storage axis exactly once");
                seen[rank] = true;
                ranks_[i] = rank;
                reversed_[i] = (value < 0);
            }
        }

        storageDims.resize(n);
        for (std::size_t i=0; i<n; i++)
            storageDims[ranks_[i]] = dims[i];

        // Strides of the storage in its own order, and of the view in its own
        std::vector<Offset> storageStrides(n), viewStrides(n);
        Offset stride = 1;
        for (std::size_t r=0; r<n; r++)
        {
            storageStrides[r] = stride;
            stride *= static_cast<Offset>(storageDims[r]);
        }
        stride = 1;
        for (std::size_t i=0; i<n; i++)
        {
            viewStrides[i] = stride;
            stride *= static_cast<Offset>(dims[i]);
        }

        strides.resize(n);
        inverseStrides.resize(n);
        for (std::size_t i=0; i<n; i++)
        {
            const int rank = ranks_[i];
            const Offset last = (dims[i] > 0 ? static_cast<Offset>(dims[i]) - 1 : 0);
            if (reversed_[i])
            {
                strides[i] = -storageStrides[rank];
                base += last * storageStrides[rank];
                inverseStrides[rank] = -viewStrides[i];
                inverseBase += last * viewStrides[i];
            }
            else
            {
                strides[i] = storageStrides[rank];
                inverseStrides[rank] = viewStrides[i];
            }

            if (rank != static_cast<int>(i) || reversed_[i])
                identity = false;
        }
    }

    Extent size () const
    {
        Extent total = 1;
        for (std::size_t i=0; i<dims.size(); i++)
            total *= dims[i];
        return total;
    }

    // Visits the view in its own order, yielding storage offsets relative to
    // base
    OffsetWalker walker () const { return OffsetWalker(dims, strides); }

    // Visits the storage in its own order, yielding view offsets relative to
    // inverseBase. This is the direction a stream is decoded in
    OffsetWalker inverseWalker () const { return OffsetWalker(storageDims, inverseStrides); }

    // The storage offset of one zero-based view linear index
    Extent storageIndex (const Extent viewIndex) const
    {
        Offset offset = base;
        Extent remainder = viewIndex;
        for (std::size_t i=0; i<dims.size(); i++)
        {
            if (dims[i] == 0)
                return 0;
            offset += static_cast<Offset>(remainder % dims[i]) * strides[i];
            remainder /= dims[i];
        }
        return static_cast<Extent>(offset);
    }
};

// A layout argument from R, where NULL or a zero-length vector means the
// identity
inline std::vector<int> layoutFrom (SEXP layout)
{
    if (Rf_isNull(layout) || Rf_xlength(layout) == 0)
        return std::vector<int>();
    const Rcpp::IntegerVector values(layout);
    return std::vector<int>(values.begin(), values.end());
}

// Copy an image's values out in view order, from anything with an integer
// subscript over its storage
template <typename Accessor, typename T>
inline void gatherView (const Accessor &source, const ViewMap &map, T *out)
{
    OffsetWalker walker = map.walker();
    gather(source, map.base, walker, out);
}

} // namespace imply

#endif
