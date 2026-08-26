#ifndef _IMPLY_ITERATOR_H_
#define _IMPLY_ITERATOR_H_

#include <cstddef>
#include <iterator>
#include <type_traits>

namespace imply {

// A random-access iterator over a strided run of values: a raw pointer with a
// fixed step between elements. Construct one at a dense buffer's own element
// stride (1) for ordinary contiguous access, or at a Raster's stride(dim) for
// a line running along one axis of it.
//
// OffsetWalker and gather() (Blocks.h) take the opposite approach on
// purpose -- copying a line into a contiguous buffer is what lets one kernel
// serve dense, packed and sparse storage alike. StrideIterator is for the
// complementary case: the storage is already dense, an algorithm wants to
// work on it in place, and a copy would be pure waste.
//
// T may be const, for a read-only iterator over the same kind of run
template <typename T>
class StrideIterator
{
public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename std::remove_cv<T>::type value_type;
    typedef std::ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;

    StrideIterator () : ptr_(nullptr), stride_(0) {}
    explicit StrideIterator (T *ptr, const difference_type stride = 1) : ptr_(ptr), stride_(stride) {}

    reference operator* () const { return *ptr_; }
    pointer operator-> () const { return ptr_; }
    reference operator[] (const difference_type n) const { return ptr_[n * stride_]; }

    StrideIterator & operator++ () { ptr_ += stride_; return *this; }
    StrideIterator operator++ (int) { StrideIterator tmp(*this); ptr_ += stride_; return tmp; }
    StrideIterator & operator-- () { ptr_ -= stride_; return *this; }
    StrideIterator operator-- (int) { StrideIterator tmp(*this); ptr_ -= stride_; return tmp; }

    StrideIterator & operator+= (const difference_type n) { ptr_ += n * stride_; return *this; }
    StrideIterator & operator-= (const difference_type n) { ptr_ -= n * stride_; return *this; }
    StrideIterator operator+ (const difference_type n) const { return StrideIterator(ptr_ + n*stride_, stride_); }
    StrideIterator operator- (const difference_type n) const { return StrideIterator(ptr_ - n*stride_, stride_); }
    friend StrideIterator operator+ (const difference_type n, const StrideIterator &it) { return it + n; }

    // Both sides are assumed to share a stride, which always holds in practice
    // since the stride never changes after construction, and subtracting
    // unrelated iterators never makes sense
    difference_type operator- (const StrideIterator &other) const { return (ptr_ - other.ptr_) / stride_; }

    bool operator== (const StrideIterator &other) const { return ptr_ == other.ptr_; }
    bool operator!= (const StrideIterator &other) const { return ptr_ != other.ptr_; }
    bool operator<  (const StrideIterator &other) const { return ptr_ < other.ptr_; }
    bool operator>  (const StrideIterator &other) const { return ptr_ > other.ptr_; }
    bool operator<= (const StrideIterator &other) const { return ptr_ <= other.ptr_; }
    bool operator>= (const StrideIterator &other) const { return ptr_ >= other.ptr_; }

private:
    T *ptr_;
    difference_type stride_;
};

} // namespace imply

#endif
