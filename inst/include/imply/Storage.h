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
enum class storageType
{
    logical, integer, real, complex
};

// Type tags carry both the C type and the R storage mode, which matters
// because logical and integer share a representation but not their notion of
// missingness or their arithmetic
struct logicalTag
{
    typedef int type;
    static constexpr storageType kind = storageType::logical;
    static constexpr int sexpType = LGLSXP;
    static bool isNA (const type x) { return x == NA_LOGICAL; }
    static type na () { return NA_LOGICAL; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const type x) { return x == 0; }
};

struct integerTag
{
    typedef int type;
    static constexpr storageType kind = storageType::integer;
    static constexpr int sexpType = INTSXP;
    static bool isNA (const type x) { return x == NA_INTEGER; }
    static type na () { return NA_INTEGER; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const type x) { return x == 0; }
};

struct realTag
{
    typedef double type;
    static constexpr storageType kind = storageType::real;
    static constexpr int sexpType = REALSXP;
    static bool isNA (const type x) { return ISNAN(x); }
    static type na () { return NA_REAL; }
    // NA is deliberately not zero, so a location holding one is kept when
    // an image is packed rather than being folded away
    static bool isZero (const type x) { return x == 0.0; }
};

struct complexTag
{
    typedef Rcomplex type;
    static constexpr storageType kind = storageType::complex;
    static constexpr int sexpType = CPLXSXP;
    static bool isNA (const type x) { return ISNAN(x.r) || ISNAN(x.i); }
    static type na () { Rcomplex z; z.r = NA_REAL; z.i = NA_REAL; return z; }
    static bool isZero (const type x) { return x.r == 0.0 && x.i == 0.0; }
};

// A non-owning typed window onto memory, which may belong to R or to us. This
// is what keeps a plain R array usable with no copy in either direction.
// Construct it from the pointer dispatchType() hands over, or from an Rcpp
// vector's begin()
template <typename T>
class view
{
protected:
    T *data_;
    Extent size_;

public:
    typedef T element;

    view () : data_(nullptr), size_(0) {}
    view (T *data, const Extent size) : data_(data), size_(size) {}

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
class buffer
{
protected:
    std::vector<T> data_;

public:
    typedef T element;

    buffer () {}
    explicit buffer (const Extent size, const T value = T()) : data_(size, value) {}

    T * data () { return data_.data(); }
    const T * data () const { return data_.data(); }
    Extent size () const { return data_.size(); }
    bool empty () const { return data_.empty(); }

    T & operator[] (const Extent n) { return data_[n]; }
    const T & operator[] (const Extent n) const { return data_[n]; }

    view<T> asView () { return view<T>(data_.data(), data_.size()); }
};

} // namespace imply

#endif
