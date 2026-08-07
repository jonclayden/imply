#ifndef _IMPLY_BLOCKS_H_
#define _IMPLY_BLOCKS_H_

#include <cstddef>
#include <vector>

#include "Raster.h"

namespace imply {

// Walks an arbitrary subset of an image's dimensions in R's order, with the
// first listed dimension moving fastest, yielding successive memory offsets.
//
// The offset is maintained incrementally by an odometer rather than being
// recomputed, and no index is materialised, so a traversal costs one addition
// per step in the common case and allocates nothing. That is what allows a
// sub-array to be gathered without first permuting the whole image the way
// base::apply() does
class offsetWalker
{
protected:
    std::vector<Extent> dims_, strides_, loc_;
    Offset offset_;
    Extent position_, total_;

public:
    offsetWalker () : offset_(0), position_(0), total_(0) {}

    offsetWalker (const std::vector<Extent> &dims, const std::vector<Extent> &strides)
        : dims_(dims), strides_(strides), loc_(dims.size(), 0), offset_(0), position_(0)
    {
        total_ = 1;
        for (std::size_t i=0; i<dims_.size(); i++)
            total_ *= dims_[i];
    }

    // The number of positions visited, which is one for an empty dimension set
    Extent size () const { return total_; }

    Offset offset () const { return offset_; }

    void reset () { seek(0); }

    // Jump straight to a position, so that a chunk of work can start part way
    // through without walking everything before it
    void seek (const Extent position)
    {
        position_ = position;
        offset_ = 0;

        Extent remainder = position;
        for (std::size_t i=0; i<dims_.size(); i++)
        {
            loc_[i] = remainder % dims_[i];
            offset_ += static_cast<Offset>(loc_[i] * strides_[i]);
            remainder /= dims_[i];
        }
    }

    // Advance to the next position, returning false once exhausted
    bool next ()
    {
        if (++position_ >= total_)
            return false;

        for (std::size_t i=0; i<dims_.size(); i++)
        {
            offset_ += static_cast<Offset>(strides_[i]);
            if (++loc_[i] < dims_[i])
                return true;

            // This dimension has wrapped, so rewind it and carry
            loc_[i] = 0;
            offset_ -= static_cast<Offset>(dims_[i] * strides_[i]);
        }

        return true;
    }
};

// Copy the values a walker visits into a contiguous buffer, starting from a
// base offset.
//
// The source is anything with an integer subscript, which is what lets one
// kernel serve dense, packed and sparse images alike: this gather is where a
// narrow type is widened and where an absent sparse location becomes a zero,
// so nothing downstream has to know which it was reading from
template <typename Accessor, typename T>
inline void gather (const Accessor &source, const Offset base, offsetWalker &walker, T *out)
{
    const Extent n = walker.size();
    walker.reset();
    for (Extent i=0; i<n; i++)
    {
        out[i] = source[static_cast<Extent>(base + walker.offset())];
        walker.next();
    }
}

} // namespace imply

#endif
