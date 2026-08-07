#ifndef _IMPLY_SPARSE_H_
#define _IMPLY_SPARSE_H_

#include <Rcpp.h>

#include <cstdint>
#include <cstring>
#include <vector>

#include "Raster.h"

namespace imply {

namespace internal {

inline int popcount64 (std::uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    // Portable fallback, in case a compiler without the builtin turns up
    x = x - ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
    return static_cast<int>((x * 0x0101010101010101ull) >> 56);
#endif
}

// The mask lives in an R raw vector, whose data is not guaranteed to be
// aligned for a 64-bit read, so words are copied out rather than cast. Every
// compiler worth the name turns this into a single load
inline std::uint64_t wordAt (const Rbyte *bytes, const std::size_t word)
{
    std::uint64_t value;
    std::memcpy(&value, bytes + word * 8, 8);
    return value;
}

} // namespace internal

// Which spatial locations of an image hold data, as one bit each.
//
// Sparsity is over locations, not over individual values: a location is either
// present, in which case the whole vector of values held there is stored, or
// absent, in which case every one of them is implicitly zero. That is what a
// brain mask actually is, and it means the packed values stay contiguous.
//
// Locating a value is O(1), via a prefix count of the bits set in every
// preceding word. tractor.base's SparseArray instead matches against a
// recomputed vector of linear indices, which is O(nnz) for every single index
// operation.
class locationMask
{
protected:
    const Rbyte *bytes;
    std::vector<Extent> blockRank;
    Extent locations_, count_;

public:
    locationMask (SEXP mask, const Extent locations)
        : bytes(RAW(mask)), locations_(locations), count_(0)
    {
        const std::size_t words = (locations + 63) / 64;
        if (static_cast<std::size_t>(Rf_xlength(mask)) < words * 8)
            Rcpp::stop("Mask is too short for an image with %d locations", double(locations));

        blockRank.resize(words + 1);
        for (std::size_t w=0; w<words; w++)
        {
            blockRank[w] = count_;
            count_ += static_cast<Extent>(internal::popcount64(internal::wordAt(bytes, w)));
        }
        blockRank[words] = count_;
    }

    Extent locations () const { return locations_; }

    // Number of locations present
    Extent count () const { return count_; }

    bool test (const Extent i) const
    {
        return (internal::wordAt(bytes, i >> 6) >> (i & 63)) & 1ull;
    }

    // Index of location i among those present, valid only when test(i)
    Extent rank (const Extent i) const
    {
        const std::uint64_t word = internal::wordAt(bytes, i >> 6);
        const std::uint64_t below = (1ull << (i & 63)) - 1ull;
        return blockRank[i >> 6] + static_cast<Extent>(internal::popcount64(word & below));
    }
};

// Reads an image whose locations are masked and whose values are packed, as
// though it were dense. The dense linear index is split into a location and a
// position within that location's values, which works because the spatial
// dimensions lead and are contiguous
template <typename T>
class sparseAccessor
{
protected:
    const locationMask &mask;
    const T *values;
    Extent locations, elements;
    T zero;

public:
    sparseAccessor (const locationMask &mask, const T *values, const Extent elements, const T zero = T())
        : mask(mask), values(values), locations(mask.locations()), elements(elements), zero(zero) {}

    T operator[] (const Extent n) const
    {
        const Extent location = n % locations;
        if (!mask.test(location))
            return zero;
        return values[mask.rank(location) * elements + (n / locations)];
    }
};

// A plain pointer wearing the same interface, so a kernel can be written once
// against either
template <typename T>
class denseAccessor
{
protected:
    const T *values;

public:
    explicit denseAccessor (const T *values) : values(values) {}
    T operator[] (const Extent n) const { return values[n]; }
};

} // namespace imply

#endif
