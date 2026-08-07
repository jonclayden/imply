#ifndef _IMPLY_STORAGE_H_
#define _IMPLY_STORAGE_H_

#include <Rcpp.h>

#include <cstddef>
#include <vector>

#include "Raster.h"

namespace imply {

// The storage types currently understood. R's four relevant modes come first;
// the narrow NIfTI-style types are added later, and are handled by converting
// once per block during the gather rather than by instantiating every kernel
// against every type
enum class StorageType
{
    logical, integer, real, complex
};

// Type tags carry both the C type and the R storage mode, which matters
// because logical and integer share a representation but not their notion of
// missingness or their arithmetic
struct LogicalTag
{
    typedef int Type;
    static constexpr StorageType kind = StorageType::logical;
    static constexpr int sexpType = LGLSXP;
    static bool isNA (const Type x) { return x == NA_LOGICAL; }
    static Type na () { return NA_LOGICAL; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const Type x) { return x == 0; }
};

struct IntegerTag
{
    typedef int Type;
    static constexpr StorageType kind = StorageType::integer;
    static constexpr int sexpType = INTSXP;
    static bool isNA (const Type x) { return x == NA_INTEGER; }
    static Type na () { return NA_INTEGER; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const Type x) { return x == 0; }
};

struct RealTag
{
    typedef double Type;
    static constexpr StorageType kind = StorageType::real;
    static constexpr int sexpType = REALSXP;
    static bool isNA (const Type x) { return ISNAN(x); }
    static Type na () { return NA_REAL; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const Type x) { return x == 0.0; }
};

struct ComplexTag
{
    typedef Rcomplex Type;
    static constexpr StorageType kind = StorageType::complex;
    static constexpr int sexpType = CPLXSXP;
    static bool isNA (const Type x) { return ISNAN(x.r) || ISNAN(x.i); }
    static Type na () { Rcomplex z; z.r = NA_REAL; z.i = NA_REAL; return z; }
    static bool isZero (const Type x) { return x.r == 0.0 && x.i == 0.0; }
};

// A non-owning typed window onto memory, which may belong to R or to us. This
// is what keeps a plain R array usable with no copy in either direction.
// Construct it from the pointer dispatchType() hands over, or from an Rcpp
// vector's begin()
template <typename T>
class View
{
protected:
    T *data_;
    Extent size_;

public:
    typedef T Element;

    View () : data_(nullptr), size_(0) {}
    View (T *data, const Extent size) : data_(data), size_(size) {}

    T * data () { return data_; }
    const T * data () const { return data_; }
    Extent size () const { return size_; }
    bool empty () const { return size_ == 0 || data_ == nullptr; }

    T & operator[] (const Extent n) { return data_[n]; }
    const T & operator[] (const Extent n) const { return data_[n]; }

    T * begin () { return data_; }
    T * end () { return data_ + size_; }
    const T * begin () const { return data_; }
    const T * end () const { return data_ + size_; }
};

// Owning storage, used only for intermediates that have no R counterpart
template <typename T>
class Buffer
{
protected:
    std::vector<T> data_;

public:
    typedef T Element;

    Buffer () {}
    explicit Buffer (const Extent size, const T value = T()) : data_(size, value) {}

    T * data () { return data_.data(); }
    const T * data () const { return data_.data(); }
    Extent size () const { return data_.size(); }
    bool empty () const { return data_.empty(); }

    T & operator[] (const Extent n) { return data_[n]; }
    const T & operator[] (const Extent n) const { return data_[n]; }

    View<T> asView () { return View<T>(data_.data(), data_.size()); }
};

} // namespace imply

#endif
