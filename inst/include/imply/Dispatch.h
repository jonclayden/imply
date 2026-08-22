#ifndef _IMPLY_DISPATCH_H_
#define _IMPLY_DISPATCH_H_

#include <R.h>
#include <Rinternals.h>

#include <type_traits>

#include "Raster.h"
#include "Storage.h"

namespace imply {

// Dispatch happens exactly once, at the R boundary, so that everything inside
// is fully typed and free of per-element branching or virtual calls

// Call fn with the type tag matching x's storage mode. fn is a generic lambda
// taking (Tag, Tag::Type *)
template <class Functor>
inline SEXP dispatchType (SEXP x, Functor &&fn)
{
    switch (TYPEOF(x))
    {
        case LGLSXP:  return fn(LogicalTag(), LOGICAL(x));
        case INTSXP:  return fn(IntegerTag(), INTEGER(x));
        case REALSXP: return fn(RealTag(), REAL(x));
        case CPLXSXP: return fn(ComplexTag(), COMPLEX(x));
        default:
        Rf_error("Unsupported storage mode '%s': imply handles logical, integer, double and complex data", Rf_type2char(TYPEOF(x)));
    }
}

// Call fn with a compile-time dimensionality where one of the common cases
// applies, and with the runtime-dimensionality Raster otherwise. fn is a
// generic lambda taking std::integral_constant<int,D>.
//
// The fixed variants exist because a dynamic Raster allocates when it builds an
// index, which in a per-voxel loop is an allocation per voxel. Instantiating
// beyond five dimensions has no practical payoff for image data
template <class Functor>
inline SEXP dispatchDims (const int nDims, Functor &&fn)
{
    switch (nDims)
    {
        case 1: return fn(std::integral_constant<int,1>());
        case 2: return fn(std::integral_constant<int,2>());
        case 3: return fn(std::integral_constant<int,3>());
        case 4: return fn(std::integral_constant<int,4>());
        case 5: return fn(std::integral_constant<int,5>());
        default: return fn(std::integral_constant<int,dynamic>());
    }
}

} // namespace imply

#endif
