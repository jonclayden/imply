#ifndef _IMPLY_RASTER_H_
#define _IMPLY_RASTER_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace imply {

typedef std::size_t Extent;
typedef std::ptrdiff_t Offset;

// Sentinel dimensionality indicating that the number of dimensions is only
// known at runtime
constexpr int dynamic = -1;

// A contiguous run of spatial locations. This is the unit of blocked traversal
// and, equivalently, of parallel work division
struct Block
{
    Extent start, length;

    Block () : start(0), length(0) {}
    Block (const Extent start, const Extent length) : start(start), length(length) {}

    Extent end () const { return start + length; }
};

namespace internal {

// Compile-time unrolled index flattening, drawn from tractor.track's Image.h.
// The recursion is on N, the number of dimensions still to fold in.
//
// NB: unlike the original, the base case multiplies by strides[0] rather than
// returning loc[0] directly. The original was only correct when the first
// dimension is contiguous, which is not true of a permuted view
template <int D, int N=D>
struct Indexer
{
    static Extent flatten (const std::array<Extent,D> &loc, const std::array<Extent,D> &strides)
    {
        return strides[N-1] * loc[N-1] + Indexer<D,N-1>::flatten(loc, strides);
    }
};

template <int D>
struct Indexer<D,1>
{
    static Extent flatten (const std::array<Extent,D> &loc, const std::array<Extent,D> &strides)
    {
        return strides[0] * loc[0];
    }
};

// Extents are held on the stack when the dimensionality is fixed, and on the
// heap only when it must be discovered at runtime. This is the substantive
// reason for having a compile-time variant at all: it removes an allocation
// per index construction from hot loops
template <int D> struct ExtentContainer { typedef std::array<Extent,D> Type; };
template <> struct ExtentContainer<dynamic> { typedef std::vector<Extent> Type; };

} // namespace internal

// A bounded n-dimensional index space with general strides. Carries no data.
//
// The dimension vector is split at `spatial`: leading dimensions index
// location, trailing dimensions index the value held at each location (a time
// series, vector, tensor, etc.). All indexing goes through the stride vector,
// so a permuted or sliced view is a stride permutation rather than a copy
template <int D = dynamic>
class Raster
{
public:
    typedef typename internal::ExtentContainer<D>::Type Index;
    static constexpr bool isFixed = (D != dynamic);

protected:
    Index dims_, strides_;
    Extent length_, spatialLength_, elementLength_;
    int spatial_;

    // Row-major in the R sense: the first index moves fastest
    void calculateStrides ()
    {
        const int n = static_cast<int>(dims_.size());
        Extent stride = 1;
        for (int i=0; i<n; i++)
        {
            strides_[i] = stride;
            stride *= dims_[i];
        }
    }

    void calculateLengths ()
    {
        const int n = static_cast<int>(dims_.size());
        if (spatial_ < 0 || spatial_ > n)
            throw std::runtime_error("Number of spatial dimensions is out of range");

        spatialLength_ = 1;
        for (int i=0; i<spatial_; i++)
            spatialLength_ *= dims_[i];

        elementLength_ = 1;
        for (int i=spatial_; i<n; i++)
            elementLength_ *= dims_[i];

        length_ = spatialLength_ * elementLength_;
    }

    void resizeIfDynamic (const std::size_t n)
    {
        if constexpr (!isFixed)
        {
            dims_.resize(n);
            strides_.resize(n);
        }
        else if (n != static_cast<std::size_t>(D))
            throw std::runtime_error("Dimension vector is not of the right dimensionality");
    }

public:
    Raster () : length_(0), spatialLength_(0), elementLength_(0), spatial_(0)
    {
        if constexpr (isFixed)
        {
            dims_.fill(0);
            strides_.fill(0);
        }
    }

    // Not explicit: we want a dimension vector to convert to a Raster freely
    template <typename Container>
    Raster (const Container &dims, const int spatial = -1)
    {
        resizeIfDynamic(dims.size());
        std::copy(dims.begin(), dims.end(), dims_.begin());
        spatial_ = (spatial < 0 ? std::min<int>(3, static_cast<int>(dims_.size())) : spatial);
        calculateStrides();
        calculateLengths();
    }

    // Construct a view with explicit strides, as produced by permutation or
    // slicing. No ownership or bounds relationship to any parent is implied
    template <typename Container>
    Raster (const Container &dims, const Container &strides, const int spatial)
    {
        if (dims.size() != strides.size())
            throw std::runtime_error("Dimension and stride vectors are of different lengths");
        resizeIfDynamic(dims.size());
        std::copy(dims.begin(), dims.end(), dims_.begin());
        std::copy(strides.begin(), strides.end(), strides_.begin());
        spatial_ = spatial;
        calculateLengths();
    }

    int nDims () const { return static_cast<int>(dims_.size()); }
    int spatial () const { return spatial_; }

    const Index & dim () const { return dims_; }
    const Index & strides () const { return strides_; }

    Extent dim (const int i) const { return dims_[i]; }
    Extent stride (const int i) const { return strides_[i]; }

    // Total number of elements, the number of spatial locations, and the number
    // of values held at each location
    Extent size () const { return length_; }
    Extent spatialSize () const { return spatialLength_; }
    Extent elementSize () const { return elementLength_; }

    bool empty () const { return length_ == 0; }

    // True if the whole Raster is densely packed in the default stride order,
    // which is what allows a kernel to fall back on plain pointer arithmetic
    bool isContiguous () const
    {
        Extent stride = 1;
        for (int i=0; i<nDims(); i++)
        {
            if (strides_[i] != stride)
                return false;
            stride *= dims_[i];
        }
        return true;
    }

    bool isContiguous (const int dim) const { return strides_[dim] == 1; }

    Extent flattenIndex (const Index &loc) const
    {
        if constexpr (isFixed)
            return internal::Indexer<D>::flatten(loc, strides_);
        else
        {
            Extent result = 0;
            for (int i=0; i<nDims(); i++)
                result += loc[i] * strides_[i];
            return result;
        }
    }

    void expandIndex (const Extent n, Index &result) const
    {
        // Only meaningful in the default stride order; a permuted view must be
        // expanded against its own dimension order
        Extent remainder = n;
        for (int i=0; i<nDims(); i++)
        {
            // A zero-extent dimension has nothing to divide by; see the same
            // guard in OffsetWalker::seek() (Blocks.h)
            if (dims_[i] == 0)
            {
                result[i] = 0;
                continue;
            }
            result[i] = remainder % dims_[i];
            remainder /= dims_[i];
        }
    }

    Index expandIndex (const Extent n) const
    {
        Index result;
        if constexpr (!isFixed)
            result.resize(dims_.size());
        expandIndex(n, result);
        return result;
    }

    // Offset of the eth value within the block of values at a single location.
    // The fast path covers any element sub-space that is itself densely packed,
    // which is every case arising from an R array
    Offset elementOffset (const Extent e) const
    {
        const int n = nDims();
        if (spatial_ >= n)
            return 0;

        bool packed = true;
        Extent stride = strides_[spatial_];
        for (int i=spatial_; i<n; i++)
        {
            if (strides_[i] != stride)
            {
                packed = false;
                break;
            }
            stride *= dims_[i];
        }

        if (packed)
            return static_cast<Offset>(e * strides_[spatial_]);

        Offset result = 0;
        Extent remainder = e;
        for (int i=spatial_; i<n; i++)
        {
            // Zero-extent dimension: nothing to divide by, see seek() in Blocks.h
            if (dims_[i] == 0)
                continue;
            result += static_cast<Offset>((remainder % dims_[i]) * strides_[i]);
            remainder /= dims_[i];
        }
        return result;
    }

    // Offset of a spatial location given its linear index among locations
    Offset spatialOffset (const Extent n) const
    {
        Offset result = 0;
        Extent remainder = n;
        for (int i=0; i<spatial_; i++)
        {
            // Zero-extent dimension: nothing to divide by, see seek() in Blocks.h
            if (dims_[i] == 0)
                continue;
            result += static_cast<Offset>((remainder % dims_[i]) * strides_[i]);
            remainder /= dims_[i];
        }
        return result;
    }

    // Number of one-dimensional lines running along `dim`, and the offset of
    // the nth of them. Lines never overlap, so they can be processed
    // independently without synchronisation
    Extent countLines (const int dim) const
    {
        Extent n = 1;
        for (int i=0; i<nDims(); i++)
        {
            if (i != dim)
                n *= dims_[i];
        }
        return n;
    }

    Offset lineOffset (const Extent n, const int dim) const
    {
        Offset result = 0;
        Extent remainder = n;
        for (int i=0; i<nDims(); i++)
        {
            if (i == dim)
                continue;
            // Zero-extent dimension: nothing to divide by, see seek() in Blocks.h
            if (dims_[i] == 0)
                continue;
            // The usual stride doesn't apply because one dimension is skipped
            result += static_cast<Offset>((remainder % dims_[i]) * strides_[i]);
            remainder /= dims_[i];
        }
        return result;
    }

    // Partition the spatial locations into runs small enough that one run's
    // worth of values stays in cache. This same partition is what gets handed
    // to worker threads
    std::vector<Block> blocks (const Extent targetElements) const
    {
        const Extent perLocation = std::max<Extent>(1, elementLength_);
        Extent size = std::max<Extent>(1, targetElements / perLocation);
        return blocksOfSize(size);
    }

    std::vector<Block> blocksOfSize (const Extent size) const
    {
        std::vector<Block> result;
        if (spatialLength_ == 0 || size == 0)
            return result;

        result.reserve((spatialLength_ + size - 1) / size);
        for (Extent start=0; start<spatialLength_; start+=size)
            result.emplace_back(start, std::min(size, spatialLength_ - start));
        return result;
    }

    // Split into at most `count` roughly equal runs. Dispatching over these
    // rather than over raw iterations is what makes a requested thread count
    // meaningful on backends that offer no width control of their own
    std::vector<Block> blocksForCount (const Extent count) const
    {
        if (count == 0 || spatialLength_ == 0)
            return std::vector<Block>();
        return blocksOfSize((spatialLength_ + count - 1) / count);
    }

    // A permuted view. Both the dimension and stride vectors are reordered, so
    // no data movement is implied
    template <typename Container>
    Raster<D> permute (const Container &order) const
    {
        if (static_cast<int>(order.size()) != nDims())
            throw std::runtime_error("Permutation is not of the right length");

        Index newDims = dims_, newStrides = strides_;
        std::vector<bool> seen(nDims(), false);
        for (int i=0; i<nDims(); i++)
        {
            const int j = static_cast<int>(order[i]);
            if (j < 0 || j >= nDims() || seen[j])
                throw std::runtime_error("Permutation is not a valid ordering");
            seen[j] = true;
            newDims[i] = dims_[j];
            newStrides[i] = strides_[j];
        }

        return Raster<D>(newDims, newStrides, spatial_);
    }
};

typedef Raster<dynamic> DynamicRaster;

template <int D>
using FixedRaster = Raster<D>;

} // namespace imply

#endif
