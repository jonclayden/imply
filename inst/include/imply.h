#ifndef _IMPLY_H_
#define _IMPLY_H_

// Public C++ interface to imply.
//
// Everything is header-only, so a package need only declare
//
//     LinkingTo: imply
//
// in its DESCRIPTION and include this file. There is nothing to link against.
//
// The pieces most likely to be wanted from outside are:
//
//   Raster<D>        an n-dimensional index space with general strides, split
//                    into spatial and value dimensions. Use FixedRaster<D>
//                    where the dimensionality is known when compiling, which
//                    keeps extents and indices on the stack, and DynamicRaster
//                    where it is not
//
//   ImageSpace       voxel-to-world geometry, with no dependency on any file
//                    format, plus the point conversions and rounding
//                    strategies that go with it
//
//   OffsetWalker     traversal of an arbitrary subset of dimensions, yielding
//                    memory offsets without materialising an index
//
//   StrideIterator   a random-access iterator over a strided run of values,
//                    for handing a dense line straight to a standard
//                    algorithm rather than gathering it into a buffer first
//
//   parallelFor      work division over libdispatch, OpenMP or neither
//
//   LocationMask     a bitset over spatial locations with O(1) rank, and the
//   SparseAccessor   accessor that reads packed values as though dense
//
//   NarrowAccessor   reads NIfTI-style narrow storage as double
//
// The accessors share one interface deliberately: a kernel written against
// it serves dense, sparse and narrowly-stored images alike.

#include "imply/Raster.h"
#include "imply/Space.h"
#include "imply/Blocks.h"
#include "imply/Iterator.h"
#include "imply/Parallel.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Sparse.h"
#include "imply/Narrow.h"
#include "imply/Sink.h"

// Linked to the package version as 100 * (major version) + (minor version). May not
// change if the API does not change, and in particular never changes with patch level
#define IMPLY_API_VERSION 2

#endif
