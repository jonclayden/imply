#ifndef _IMPLY_SINK_H_
#define _IMPLY_SINK_H_

#include <Rcpp.h>

#include <memory>

namespace imply {

// True if a value carries anything beyond its raw data. The fast path copies
// values into a flat vector, which would discard such attributes, so anything
// carrying them is diverted to a list and shaped in R instead.
//
// Rcpp reports every attribute, where enumerating them through the public C
// API would mean naming each one and quietly missing the rest. An object with
// no attributes returns an empty vector without allocating
inline bool hasAttributes (SEXP x)
{
    return !Rcpp::RObject(x).attributeNames().empty();
}

// Results are written through a sink rather than accumulated in a list and
// simplified afterwards. When every call returns the same shape and type --
// which is the case for any reduction -- the output vector is allocated once
// up front and filled in place, so nothing proportional to the number of calls
// is held as separate R objects.
//
// The abstraction also leaves room for a file-backed sink, for results too
// large to hold in memory, without any kernel having to change.
class Sink
{
public:
    virtual ~Sink () {}

    // Returns false if this sink cannot represent the value, in which case the
    // caller falls back to a list
    virtual bool write (const R_xlen_t index, SEXP value) = 0;

    virtual SEXP finish () = 0;
};

// The fast path: a preallocated atomic vector holding results of uniform type
// and length
class VectorSink : public Sink
{
public:
    virtual R_xlen_t length () const = 0;

    // Unpack what has been written so far, so that a fallback list can take
    // over mid-run without the completed calls being repeated
    virtual SEXP element (const R_xlen_t index) const = 0;
};

// Templating on the storage type means the element copies are ordinary typed
// iterator copies, rather than a switch over R's accessor macros repeated at
// every use
template <int RTYPE>
class TypedSink : public VectorSink
{
protected:
    Rcpp::Vector<RTYPE> values;
    R_xlen_t elementLength;

public:
    TypedSink (const R_xlen_t elementLength, const R_xlen_t count)
        : values(elementLength * count), elementLength(elementLength) {}

    R_xlen_t length () const { return elementLength; }

    bool write (const R_xlen_t index, SEXP value)
    {
        if (TYPEOF(value) != RTYPE || Rf_xlength(value) != elementLength || hasAttributes(value))
            return false;

        // Types match, so this wraps the existing vector rather than copying it
        const Rcpp::Vector<RTYPE> source(value);
        std::copy(source.begin(), source.end(), values.begin() + index * elementLength);
        return true;
    }

    SEXP element (const R_xlen_t index) const
    {
        Rcpp::Vector<RTYPE> result(elementLength);
        std::copy(values.begin() + index * elementLength,
                  values.begin() + (index + 1) * elementLength,
                  result.begin());
        return result;
    }

    SEXP finish () { return values; }
};

// The one place a runtime storage type has to become a compile-time one
inline std::unique_ptr<VectorSink> makeVectorSink (const SEXPTYPE type, const R_xlen_t elementLength,
                                                   const R_xlen_t count)
{
    switch (type)
    {
        case LGLSXP:  return std::unique_ptr<VectorSink>(new TypedSink<LGLSXP>(elementLength, count));
        case INTSXP:  return std::unique_ptr<VectorSink>(new TypedSink<INTSXP>(elementLength, count));
        case REALSXP: return std::unique_ptr<VectorSink>(new TypedSink<REALSXP>(elementLength, count));
        case CPLXSXP: return std::unique_ptr<VectorSink>(new TypedSink<CPLXSXP>(elementLength, count));
        default:      return std::unique_ptr<VectorSink>();
    }
}

// The general path, for results that vary in type, length or structure. This
// is what base::apply() always does, so falling back to it costs nothing
// relative to the status quo
class ListSink : public Sink
{
protected:
    Rcpp::List values;

public:
    explicit ListSink (const R_xlen_t count) : values(count) {}

    bool write (const R_xlen_t index, SEXP value)
    {
        values[index] = value;
        return true;
    }

    SEXP finish () { return values; }
};

} // namespace imply

#endif
