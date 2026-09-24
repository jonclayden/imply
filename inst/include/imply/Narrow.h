#ifndef _IMPLY_NARROW_H_
#define _IMPLY_NARROW_H_

#include <Rcpp.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

#include "Raster.h"

namespace imply {

// R has no single-precision type, so a large image may pay twice the memory it
// needs and, since these passes are bandwidth-bound, roughly twice the time.
// Raw MRI is commonly 16-bit, which is four times narrower than double.
//
// Values are held in an R raw vector and reinterpreted, with an optional
// affine scaling: value = stored * slope + intercept. That is the NIfTI model,
// and it lets an integer type carry a range it could not otherwise hold.
//
// A raw vector is used rather than an external pointer so that the data is
// garbage-collected, serialises, and survives a save/load like any other R
// object.
enum class NarrowType { int8, uint8, int16, uint16, int32, float32 };

inline std::string narrowTypeName (const NarrowType type)
{
    switch (type)
    {
        case NarrowType::int8:    return "int8";
        case NarrowType::uint8:   return "uint8";
        case NarrowType::int16:   return "int16";
        case NarrowType::uint16:  return "uint16";
        case NarrowType::int32:   return "int32";
        case NarrowType::float32: return "float32";
    }
    return "unknown";
}

inline NarrowType narrowTypeFromName (const std::string &name)
{
    if (name == "int8")    return NarrowType::int8;
    if (name == "uint8")   return NarrowType::uint8;
    if (name == "int16")   return NarrowType::int16;
    if (name == "uint16")  return NarrowType::uint16;
    if (name == "int32")   return NarrowType::int32;
    if (name == "float32") return NarrowType::float32;
    Rcpp::stop("Unknown storage type \"%s\"; expected int8, uint8, int16, uint16, int32 or float32", name);
}

inline std::size_t narrowTypeSize (const NarrowType type)
{
    switch (type)
    {
        case NarrowType::int8:
        case NarrowType::uint8:   return 1;
        case NarrowType::int16:
        case NarrowType::uint16:  return 2;
        case NarrowType::int32:
        case NarrowType::float32: return 4;
    }
    return 0;
}

inline bool narrowTypeIsInteger (const NarrowType type)
{
    return type != NarrowType::float32;
}

// The representable range, as doubles, used when deciding whether a scaling is
// needed and when clamping on the way in
inline void narrowTypeRange (const NarrowType type, double &low, double &high)
{
    switch (type)
    {
        case NarrowType::int8:    low = -128.0;        high = 127.0;         break;
        case NarrowType::uint8:   low = 0.0;           high = 255.0;         break;
        case NarrowType::int16:   low = -32768.0;      high = 32767.0;       break;
        case NarrowType::uint16:  low = 0.0;           high = 65535.0;       break;
        case NarrowType::int32:   low = -2147483648.0; high = 2147483647.0;  break;
        case NarrowType::float32: low = -3.4028234663852886e38;
                                  high = 3.4028234663852886e38;             break;
    }
}

namespace internal {

// Raw vector data carries no alignment guarantee, so stored values are copied
// in and out rather than dereferenced through a cast
template <typename T>
inline T readAs (const Rbyte *bytes, const Extent index)
{
    T value;
    std::memcpy(&value, bytes + index * sizeof(T), sizeof(T));
    return value;
}

template <typename T>
inline void writeAs (Rbyte *bytes, const Extent index, const T value)
{
    std::memcpy(bytes + index * sizeof(T), &value, sizeof(T));
}

} // namespace internal

// Reads narrow storage as though it were double, applying the scaling. This is
// the accessor the block gather uses, so a kernel never sees a narrow type at
// all and the number of template instantiations does not multiply by the
// number of storage types
template <typename Stored>
class NarrowAccessor
{
protected:
    const Rbyte *bytes;
    double slope, intercept;
    bool scaled;

public:
    NarrowAccessor (const Rbyte *bytes, const double slope, const double intercept)
        : bytes(bytes), slope(slope), intercept(intercept),
          scaled(slope != 1.0 || intercept != 0.0) {}

    double operator[] (const Extent n) const
    {
        const Stored stored = internal::readAs<Stored>(bytes, n);

        // A float NaN stands for a missing value; R's NA has a particular NaN
        // payload that does not survive the narrowing, so it is restored here
        if (std::numeric_limits<Stored>::has_quiet_NaN && stored != stored)
            return NA_REAL;

        const double value = static_cast<double>(stored);
        return scaled ? value * slope + intercept : value;
    }
};

// Run fn with the C type matching a narrow storage type, so a loop over
// narrow data is written once rather than six times
template <class Functor>
inline SEXP dispatchNarrowType (const NarrowType type, Functor &&fn)
{
    switch (type)
    {
        case NarrowType::int8:    return fn(std::int8_t());
        case NarrowType::uint8:   return fn(std::uint8_t());
        case NarrowType::int16:   return fn(std::int16_t());
        case NarrowType::uint16:  return fn(std::uint16_t());
        case NarrowType::int32:   return fn(std::int32_t());
        case NarrowType::float32: return fn(float());
    }
    Rcpp::stop("Unhandled storage type");
}

// Choose a scaling that maps the data onto the whole of an integer type's
// range, but only when the data does not already fit. Modelled on RNifti's
// NiftiImageData::calibrateFrom
struct Calibration
{
    double slope, intercept;
};

inline Calibration calibrateFor (const NarrowType type, const double dataMin, const double dataMax,
                                 const bool integral)
{
    Calibration result { 1.0, 0.0 };

    if (!narrowTypeIsInteger(type))
        return result;

    double low, high;
    narrowTypeRange(type, low, high);

    // R reserves the most negative int for NA, so a stored int32 of that value
    // is read as missing by anything that holds int32 values as R integers,
    // RNifti included. A scaling never needs it, so it is kept out of use
    if (type == NarrowType::int32)
        low += 1.0;

    // Whole numbers already inside the range need no scaling, and leaving it
    // alone keeps the stored values readable and exact
    if (integral && dataMin >= low && dataMax <= high)
        return result;

    if (dataMax == dataMin)
    {
        result.intercept = dataMin;
        return result;
    }

    result.slope = (dataMax - dataMin) / (high - low);
    result.intercept = dataMin - result.slope * low;
    return result;
}

} // namespace imply

#endif
