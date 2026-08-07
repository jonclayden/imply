#include <Rcpp.h>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Sparse.h"

using namespace imply;

namespace {

std::size_t maskWords (const Extent locations) { return (locations + 63) / 64; }

Rcpp::RawVector emptyMask (const Extent locations)
{
    Rcpp::RawVector mask(static_cast<R_xlen_t>(maskWords(locations) * 8));
    std::fill(mask.begin(), mask.end(), Rbyte(0));
    return mask;
}

inline void setBit (Rbyte *bytes, const Extent i)
{
    bytes[i >> 3] |= static_cast<Rbyte>(1u << (i & 7));
}

inline bool getBit (const Rbyte *bytes, const Extent i)
{
    return (bytes[i >> 3] >> (i & 7)) & 1u;
}

// Split the dense dimensions into a count of locations and a count of values
// held at each. The spatial dimensions lead and are contiguous, which is what
// lets a dense linear index be decomposed arithmetically
struct shape
{
    Extent locations, elements;
};

shape shapeOf (const std::vector<Extent> &dims, const int spatial)
{
    shape result;
    result.locations = 1;
    for (int i=0; i<spatial; i++)
        result.locations *= dims[i];
    result.elements = 1;
    for (std::size_t i=spatial; i<dims.size(); i++)
        result.elements *= dims[i];
    return result;
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List denseToSparse (Rcpp::RObject x, Rcpp::Nullable<Rcpp::IntegerVector> spatial = R_NilValue)
{
    const rasterSpec spec = specOf(x, spatial);
    const shape s = shapeOf(spec.dims, spec.spatial);

    Rcpp::List result;

    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;

        // A location is kept if any of the values held there is not exactly
        // zero. NA is not zero, so a location holding one survives packing
        Rcpp::RawVector mask = emptyMask(s.locations);
        Rbyte * const bytes = mask.begin();
        Extent count = 0;

        for (Extent i=0; i<s.locations; i++)
        {
            bool occupied = false;
            for (Extent e=0; e<s.elements; e++)
            {
                if (!Tag::isZero(data[i + e * s.locations]))
                {
                    occupied = true;
                    break;
                }
            }
            if (occupied)
            {
                setBit(bytes, i);
                count++;
            }
        }

        Rcpp::Vector<Tag::sexpType> values(static_cast<R_xlen_t>(count * s.elements));
        Extent position = 0;
        for (Extent i=0; i<s.locations; i++)
        {
            if (!getBit(bytes, i))
                continue;
            for (Extent e=0; e<s.elements; e++)
                values[position * s.elements + e] = data[i + e * s.locations];
            position++;
        }

        result = Rcpp::List::create(Rcpp::Named("mask") = mask,
                                    Rcpp::Named("values") = values,
                                    Rcpp::Named("count") = static_cast<double>(count));
        return R_NilValue;
    });

    return result;
}

// [[Rcpp::export]]
SEXP sparseToDense (Rcpp::RawVector mask, Rcpp::RObject values, Rcpp::IntegerVector dim, int spatial)
{
    std::vector<Extent> dims(dim.begin(), dim.end());
    const shape s = shapeOf(dims, spatial);
    const locationMask bits(mask, s.locations);

    return dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;

        Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(s.locations * s.elements));
        std::fill(result.begin(), result.end(), typename Tag::type());

        for (Extent i=0; i<s.locations; i++)
        {
            if (!bits.test(i))
                continue;
            const Extent position = bits.rank(i);
            for (Extent e=0; e<s.elements; e++)
                result[i + e * s.locations] = packed[position * s.elements + e];
        }

        result.attr("dim") = dim;
        return result;
    });
}

// Values at one-based dense linear indices, without materialising the image
// [[Rcpp::export]]
SEXP sparseElements (Rcpp::RawVector mask, Rcpp::RObject values, Rcpp::IntegerVector dim, int spatial,
                     Rcpp::NumericVector indices)
{
    std::vector<Extent> dims(dim.begin(), dim.end());
    const shape s = shapeOf(dims, spatial);
    const locationMask bits(mask, s.locations);
    const Extent total = s.locations * s.elements;

    return dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;

        const sparseAccessor<typename Tag::type> accessor(bits, packed, s.elements);
        Rcpp::Vector<Tag::sexpType> result(indices.size());

        for (R_xlen_t k=0; k<indices.size(); k++)
        {
            const double index = indices[k];
            if (Rcpp::NumericVector::is_na(index) || index < 1 || index > static_cast<double>(total))
                Rcpp::stop("Index %d is out of range", double(k + 1));
            result[k] = accessor[static_cast<Extent>(index) - 1];
        }

        return result;
    });
}

// [[Rcpp::export]]
double maskCount (Rcpp::RawVector mask, double locations)
{
    const locationMask bits(mask, static_cast<Extent>(locations));
    return static_cast<double>(bits.count());
}

// [[Rcpp::export]]
Rcpp::LogicalVector maskToLogical (Rcpp::RawVector mask, double locations)
{
    const Extent n = static_cast<Extent>(locations);
    Rcpp::LogicalVector result(static_cast<R_xlen_t>(n));
    const Rbyte * const bytes = mask.begin();
    for (Extent i=0; i<n; i++)
        result[i] = getBit(bytes, i);
    return result;
}

// [[Rcpp::export]]
Rcpp::RawVector maskFromLogical (Rcpp::LogicalVector present)
{
    const Extent n = static_cast<Extent>(present.size());
    Rcpp::RawVector mask = emptyMask(n);
    Rbyte * const bytes = mask.begin();
    for (Extent i=0; i<n; i++)
    {
        if (present[i] == NA_LOGICAL)
            Rcpp::stop("Mask must not contain missing values");
        if (present[i])
            setBit(bytes, i);
    }
    return mask;
}

// Combine two masks. Union is what an operation like addition needs, since a
// location present in either operand may hold a non-zero result; intersection
// is what multiplication needs, since anything else is certainly zero
// [[Rcpp::export]]
Rcpp::RawVector maskCombine (Rcpp::RawVector first, Rcpp::RawVector second, std::string how)
{
    if (first.size() != second.size())
        Rcpp::stop("Masks describe images of different sizes");

    Rcpp::RawVector result(first.size());
    if (how == "union")
    {
        for (R_xlen_t i=0; i<first.size(); i++)
            result[i] = first[i] | second[i];
    }
    else if (how == "intersection")
    {
        for (R_xlen_t i=0; i<first.size(); i++)
            result[i] = first[i] & second[i];
    }
    else
        Rcpp::stop("Masks can be combined by \"union\" or \"intersection\", not \"%s\"", how);

    return result;
}

// Drop locations whose values have all become zero, restoring the invariant
// that the mask holds exactly the non-zero locations. An operation such as
// multiplying by zero can break it, and leaving it broken would make
// sparseness() report a figure that is no longer true
// [[Rcpp::export]]
Rcpp::List tightenMask (Rcpp::RawVector mask, Rcpp::RObject values, double locations, double elements)
{
    const Extent n = static_cast<Extent>(locations);
    const Extent e = static_cast<Extent>(elements);
    const locationMask bits(mask, n);

    Rcpp::List result;

    dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;

        Rcpp::RawVector kept = emptyMask(n);
        Rbyte * const bytes = kept.begin();
        Extent count = 0;

        for (Extent i=0; i<n; i++)
        {
            if (!bits.test(i))
                continue;

            const Extent position = bits.rank(i);
            bool occupied = false;
            for (Extent k=0; k<e; k++)
            {
                if (!Tag::isZero(packed[position * e + k]))
                {
                    occupied = true;
                    break;
                }
            }

            if (occupied)
            {
                setBit(bytes, i);
                count++;
            }
        }

        Rcpp::Vector<Tag::sexpType> tightened(static_cast<R_xlen_t>(count * e));
        Extent target = 0;
        for (Extent i=0; i<n; i++)
        {
            if (!getBit(bytes, i))
                continue;
            const Extent position = bits.rank(i);
            for (Extent k=0; k<e; k++)
                tightened[target * e + k] = packed[position * e + k];
            target++;
        }

        result = Rcpp::List::create(Rcpp::Named("mask") = kept,
                                    Rcpp::Named("values") = tightened,
                                    Rcpp::Named("count") = static_cast<double>(count));
        return R_NilValue;
    });

    return result;
}

// Rewrite packed values so they line up with a different mask, filling zeros
// where the old mask held nothing and dropping what the new one excludes
// [[Rcpp::export]]
SEXP repackValues (Rcpp::RawVector oldMask, Rcpp::RObject values, Rcpp::RawVector newMask,
                   double locations, double elements)
{
    const Extent n = static_cast<Extent>(locations);
    const Extent e = static_cast<Extent>(elements);
    const locationMask before(oldMask, n);
    const locationMask after(newMask, n);

    return dispatchType(values, [&](auto tag, auto *packed) -> SEXP {
        typedef decltype(tag) Tag;

        Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(after.count() * e));
        std::fill(result.begin(), result.end(), typename Tag::type());

        for (Extent i=0; i<n; i++)
        {
            if (!after.test(i) || !before.test(i))
                continue;
            const Extent from = before.rank(i);
            const Extent to = after.rank(i);
            for (Extent k=0; k<e; k++)
                result[to * e + k] = packed[from * e + k];
        }

        return result;
    });
}
