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
//   raster<D>        an n-dimensional index space with general strides, split
//                    into spatial and value dimensions. Use fixedRaster<D>
//                    where the dimensionality is known when compiling, which
//                    keeps extents and indices on the stack, and dynamicRaster
//                    where it is not
//
//   imageSpace       voxel-to-world geometry, with no dependency on any file
//                    format, plus the point conversions and rounding
//                    strategies that go with it
//
//   offsetWalker     traversal of an arbitrary subset of dimensions, yielding
//                    memory offsets without materialising an index
//
//   parallelFor      work division over libdispatch, OpenMP or neither
//
//   locationMask     a bitset over spatial locations with O(1) rank, and the
//   sparseAccessor   accessor that reads packed values as though dense
//
//   narrowAccessor   reads NIfTI-style narrow storage as double
//
// The accessors share one interface deliberately: a kernel written against
// it serves dense, sparse and narrowly-stored images alike.

#include "imply/Raster.h"
#include "imply/Space.h"
#include "imply/Blocks.h"
#include "imply/Parallel.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/RImage.h"
#include "imply/Sparse.h"
#include "imply/Narrow.h"
#include "imply/Sink.h"

#endif
