#include <Rcpp.h>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Narrow.h"

using namespace imply;

// Range and integrality of the data, which is what decides whether a scaling
// is needed. Missing values are ignored here and reported separately
// [[Rcpp::export]]
Rcpp::List valueRange (Rcpp::RObject x)
{
    double low = R_PosInf, high = R_NegInf;
    bool integral = true, missing = false, any = false;

    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;

        if constexpr (Tag::kind == StorageType::complex)
            Rcpp::stop("Complex data cannot be stored in a narrow type");
        else
        {
            const R_xlen_t n = Rf_xlength(x);
            for (R_xlen_t i=0; i<n; i++)
            {
                const typename Tag::Type value = data[i];
                if (Tag::isNA(value))
                {
                    missing = true;
                    continue;
                }

                const double asDouble = static_cast<double>(value);
                if (!R_FINITE(asDouble))
                {
                    missing = true;
                    continue;
                }

                any = true;
                if (asDouble < low) low = asDouble;
                if (asDouble > high) high = asDouble;
                if (integral && asDouble != std::floor(asDouble))
                    integral = false;
            }
        }
        return R_NilValue;
    });

    if (!any)
    {
        low = 0.0;
        high = 0.0;
    }

    return Rcpp::List::create(Rcpp::Named("low") = low,
                              Rcpp::Named("high") = high,
                              Rcpp::Named("integral") = integral,
                              Rcpp::Named("missing") = missing);
}

// [[Rcpp::export]]
Rcpp::List calibrateStorage (std::string type, double low, double high, bool integral)
{
    const Calibration result = calibrateFor(narrowTypeFromName(type), low, high, integral);
    return Rcpp::List::create(Rcpp::Named("slope") = result.slope,
                              Rcpp::Named("intercept") = result.intercept);
}

// Convert native R data into narrow storage, held as a raw vector
// [[Rcpp::export]]
Rcpp::RawVector packNarrow (Rcpp::RObject x, std::string type, double slope = 1, double intercept = 0)
{
    const NarrowType target = narrowTypeFromName(type);
    const R_xlen_t n = Rf_xlength(x);

    if (slope == 0)
        Rcpp::stop("Storage slope must not be zero");

    Rcpp::RawVector result(n * static_cast<R_xlen_t>(narrowTypeSize(target)));
    Rbyte * const bytes = result.begin();

    double low, high;
    narrowTypeRange(target, low, high);
    const bool isInteger = narrowTypeIsInteger(target);

    dispatchType(x, [&](auto tag, auto *data) -> SEXP {
        typedef decltype(tag) Tag;

        if constexpr (Tag::kind == StorageType::complex)
            Rcpp::stop("Complex data cannot be stored in a narrow type");
        else
        {
            dispatchNarrowType(target, [&](auto stored) -> SEXP {
                typedef decltype(stored) Stored;

                for (R_xlen_t i=0; i<n; i++)
                {
                    const typename Tag::Type raw = data[i];
                    double value;

                    if (Tag::isNA(raw))
                    {
                        // Only a floating point type can carry missingness;
                        // packing to an integer type would silently lose it,
                        // so that combination is refused before we get here
                        internal::writeAs<Stored>(bytes, i,
                            static_cast<Stored>(std::numeric_limits<double>::quiet_NaN()));
                        continue;
                    }

                    value = (static_cast<double>(raw) - intercept) / slope;

                    if (isInteger)
                    {
                        value = std::nearbyint(value);
                        if (value < low) value = low;
                        if (value > high) value = high;
                    }

                    internal::writeAs<Stored>(bytes, i, static_cast<Stored>(value));
                }
                return R_NilValue;
            });
        }
        return R_NilValue;
    });

    return result;
}

// Materialise narrow storage back into a double vector
// [[Rcpp::export]]
Rcpp::NumericVector unpackNarrow (Rcpp::RawVector packed, std::string type, double count,
                                  double slope = 1, double intercept = 0)
{
    const NarrowType target = narrowTypeFromName(type);
    const R_xlen_t n = static_cast<R_xlen_t>(count);

    if (packed.size() < n * static_cast<R_xlen_t>(narrowTypeSize(target)))
        Rcpp::stop("Packed data is too short for %d values of type %s", double(n), type);

    Rcpp::NumericVector result(n);
    const Rbyte * const bytes = packed.begin();

    dispatchNarrowType(target, [&](auto stored) -> SEXP {
        typedef decltype(stored) Stored;
        const NarrowAccessor<Stored> accessor(bytes, slope, intercept);
        for (R_xlen_t i=0; i<n; i++)
            result[i] = accessor[static_cast<Extent>(i)];
        return R_NilValue;
    });

    return result;
}

// Values at one-based linear indices, without materialising the whole image
// [[Rcpp::export]]
Rcpp::NumericVector narrowElements (Rcpp::RawVector packed, std::string type, double count,
                                    Rcpp::NumericVector indices, double slope = 1, double intercept = 0)
{
    const NarrowType target = narrowTypeFromName(type);
    const Extent total = static_cast<Extent>(count);
    Rcpp::NumericVector result(indices.size());
    const Rbyte * const bytes = packed.begin();

    dispatchNarrowType(target, [&](auto stored) -> SEXP {
        typedef decltype(stored) Stored;
        const NarrowAccessor<Stored> accessor(bytes, slope, intercept);

        for (R_xlen_t k=0; k<indices.size(); k++)
        {
            const double index = indices[k];
            if (Rcpp::NumericVector::is_na(index) || index < 1 || index > static_cast<double>(total))
                Rcpp::stop("Index %d is out of range", double(k + 1));
            result[k] = accessor[static_cast<Extent>(index) - 1];
        }
        return R_NilValue;
    });

    return result;
}

// Summaries are computed in double regardless of how narrowly the values are
// stored, so accumulated error does not depend on the storage type
// [[Rcpp::export]]
Rcpp::List narrowSummary (Rcpp::RawVector packed, std::string type, double count,
                          double slope = 1, double intercept = 0, bool naRm = false)
{
    const NarrowType target = narrowTypeFromName(type);
    const Extent n = static_cast<Extent>(count);

    double total = 0.0, low = R_PosInf, high = R_NegInf;
    bool missing = false;
    Extent used = 0;

    const Rbyte * const bytes = packed.begin();

    dispatchNarrowType(target, [&](auto stored) -> SEXP {
        typedef decltype(stored) Stored;
        const NarrowAccessor<Stored> accessor(bytes, slope, intercept);

        for (Extent i=0; i<n; i++)
        {
            const double value = accessor[i];
            if (ISNAN(value))
            {
                missing = true;
                continue;
            }
            total += value;
            if (value < low) low = value;
            if (value > high) high = value;
            used++;
        }
        return R_NilValue;
    });

    const bool spoiled = missing && !naRm;

    return Rcpp::List::create(
        Rcpp::Named("sum") = (spoiled ? NA_REAL : total),
        Rcpp::Named("min") = (spoiled || used == 0 ? NA_REAL : low),
        Rcpp::Named("max") = (spoiled || used == 0 ? NA_REAL : high),
        Rcpp::Named("mean") = (spoiled || used == 0 ? NA_REAL : total / static_cast<double>(used)),
        Rcpp::Named("count") = static_cast<double>(used),
        Rcpp::Named("missing") = missing);
}
