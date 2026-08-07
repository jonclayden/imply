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

// R has no single-precision type, so a large image pays twice the memory it
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
enum class narrowType { int8, uint8, int16, uint16, int32, float32 };

inline std::string narrowTypeName (const narrowType type)
{
    switch (type)
    {
        case narrowType::int8:    return "int8";
        case narrowType::uint8:   return "uint8";
        case narrowType::int16:   return "int16";
        case narrowType::uint16:  return "uint16";
        case narrowType::int32:   return "int32";
        case narrowType::float32: return "float32";
    }
    return "unknown";
}

inline narrowType narrowTypeFromName (const std::string &name)
{
    if (name == "int8")    return narrowType::int8;
    if (name == "uint8")   return narrowType::uint8;
    if (name == "int16")   return narrowType::int16;
    if (name == "uint16")  return narrowType::uint16;
    if (name == "int32")   return narrowType::int32;
    if (name == "float32") return narrowType::float32;
    Rcpp::stop("Unknown storage type \"%s\"; expected int8, uint8, int16, uint16, int32 or float32", name);
}

inline std::size_t narrowTypeSize (const narrowType type)
{
    switch (type)
    {
        case narrowType::int8:
        case narrowType::uint8:   return 1;
        case narrowType::int16:
        case narrowType::uint16:  return 2;
        case narrowType::int32:
        case narrowType::float32: return 4;
    }
    return 0;
}

inline bool narrowTypeIsInteger (const narrowType type)
{
    return type != narrowType::float32;
}

// The representable range, as doubles, used when deciding whether a scaling is
// needed and when clamping on the way in
inline void narrowTypeRange (const narrowType type, double &low, double &high)
{
    switch (type)
    {
        case narrowType::int8:    low = -128.0;        high = 127.0;         break;
        case narrowType::uint8:   low = 0.0;           high = 255.0;         break;
        case narrowType::int16:   low = -32768.0;      high = 32767.0;       break;
        case narrowType::uint16:  low = 0.0;           high = 65535.0;       break;
        case narrowType::int32:   low = -2147483648.0; high = 2147483647.0;  break;
        case narrowType::float32: low = -3.4028234663852886e38;
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
class narrowAccessor
{
protected:
    const Rbyte *bytes;
    double slope, intercept;
    bool scaled;

public:
    narrowAccessor (const Rbyte *bytes, const double slope, const double intercept)
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
inline SEXP dispatchNarrowType (const narrowType type, Functor &&fn)
{
    switch (type)
    {
        case narrowType::int8:    return fn(std::int8_t());
        case narrowType::uint8:   return fn(std::uint8_t());
        case narrowType::int16:   return fn(std::int16_t());
        case narrowType::uint16:  return fn(std::uint16_t());
        case narrowType::int32:   return fn(std::int32_t());
        case narrowType::float32: return fn(float());
    }
    Rcpp::stop("Unhandled storage type");
}

// Choose a scaling that maps the data onto the whole of an integer type's
// range, but only when the data does not already fit. Modelled on RNifti's
// NiftiImageData::calibrateFrom
struct calibration
{
    double slope, intercept;
};

inline calibration calibrateFor (const narrowType type, const double dataMin, const double dataMax,
                                 const bool integral)
{
    calibration result { 1.0, 0.0 };

    if (!narrowTypeIsInteger(type))
        return result;

    double low, high;
    narrowTypeRange(type, low, high);

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
