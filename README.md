

[![CRAN version](http://www.r-pkg.org/badges/version/imply)](https://cran.r-project.org/package=imply) [![CI](https://github.com/jonclayden/imply/actions/workflows/ci.yaml/badge.svg)](https://github.com/jonclayden/imply/actions/workflows/ci.yaml) [![status](https://tinyverse.netlify.app/badge/imply)](https://tinyverse.netlify.app)

# imply: Efficient Generalised Images for R

The `imply` package provides data structures for *generalised images* – two- or three-dimensional rasters that may hold a vector or time series at each location rather than a single intensity — together with infrastructure for applying functions to them in a memory-efficient and optionally parallelised fashion. The primary use case is medical imaging, but nothing in the package assumes this explicitly. The name is a portmanteau of "image" and "apply", and follows in the footsteps of influential packages like [`plyr`](https://cran.r-project.org/package=plyr).

This is a general infrastructure package that has no dependency on any specific file format, although there is support for image geometry characteristics such as pixel/voxel sizes and transforms to represent a real-world spatial embedding. Reading, writing and interoperating with other image classes is largely out of scope,
and left to packages that sit above this one.

There is also no viewer, as conventions vary by application. [`RNifti::view()`](https://cran.r-project.org/package=RNifti) is one option for medical images; the standard `image()`, or [`mmand::display()`](https://cran.r-project.org/package=mmand), can be used for typical 2D bitmaps.

Key features include

- modest dependencies, including [`Rcpp`](https://cran.r-project.org/package=Rcpp) and [`RcppArray`](https://cran.r-project.org/package=RcppArray) for the C++ back-end and [`S7`](https://cran.r-project.org/package=S7) for modern R class types;
- three [storage representations](#image-types) — an ordinary R array, a masked sparse form, and a narrow packed form — served by **one** apply engine, so a kernel written
  once runs unchanged whichever way the data happens to be held;
- identical results between the package's [`imapply()` function](#applying-and-reducing-functions) and `base::apply()`, but with significantly less peak memory usage, which is particularly important for large images;
- easy [parallelism](#parallelism) using the standard `parallel` package, libdispatch (a.k.a. Grand Central Dispatch) or OpenMP, plus progress reporting; and
- a [C++ API](#c-api) exposing its core engine to other packages.

The package is in an early stage of development, and its interface may change in subsequent releases. The latest version of the package can be installed from GitHub using the `remotes` package:


``` r
## install.packages("remotes")
remotes::install_github("jonclayden/imply")
```

## Image types

`imply` supports three storage representations:

| | class | storage | good for |
|---|---|---|---|
| dense | `denseImage` | an ordinary R array, any storage mode | general use; anything not covered below |
| sparse | `sparseImage` | a bit mask over spatial locations, plus values packed in location order | images that are mostly zero outside a mask, such as data masked to a region of interest |
| packed | `packedImage` | a raw vector reinterpreted as a narrower integer or floating-point type than R usually supports | large images where the values are not all needed as `double` at once |

The back-end allows functions applied to images to target all of these image forms efficiently without specialisation for each type.

A `denseImage` wraps an array with some geometry attached: how many of the leading dimensions are spatial, the size of a voxel, and the 4x4 affine matrix mapping voxel to world coordinates. The latter can be ignored where it isn't important.


``` r
library(imply)

data <- array(rnorm(4 * 5 * 6 * 10), dim = c(4, 5, 6, 10))
image <- denseImage(data, voxelSize = c(2, 2, 2.5), spatial = 3)
image
## Dense image: 4 x 5 x 6 x 10 (double)
##   Spatial dimensions : 4 x 5 x 6
##   Voxel size         : 2 x 2 x 2.5 (unit unknown)
##   Values per location: 10
```

Here the image is 4 x 5 x 6 spatially, with a ten-point time series at each location — `spatial = 3` says so explicitly, though it would have been inferred anyway, since `spatial()` defaults to the first three dimensions or however many there are, whichever is fewer.


``` r
dim(image)
## [1]  4  5  6 10
spatial(image)
## [1] 3
voxelSize(image)
## [1] 2.0 2.0 2.5
worldTransform(image)
##      [,1] [,2] [,3] [,4]
## [1,]    2    0  0.0    0
## [2,]    0    2  0.0    0
## [3,]    0    0  2.5    0
## [4,]    0    0  0.0    1
```

Voxel size and world placement are stored, and can be set, independently of one another, so replacing one never has to touch the other — see `?geometry`. `fromVoxel()` and `toVoxel()` convert points between voxel and world space using the composed transform, using the usual R one-based coordinates for voxels.


``` r
fromVoxel(c(1, 1, 1), image)
##      [,1] [,2] [,3]
## [1,]    0    0    0
toVoxel(fromVoxel(c(1, 1, 1), image), image)
##      [,1] [,2] [,3]
## [1,]    1    1    1
```

Because a `denseImage` is a plain array underneath, ordinary R operations — indexing, arithmetic, comparison — work on it directly and return a plain array, since an arbitrary index or elementwise result has no well-defined geometry of its own.

### Converting between representations

`asSparse()`, `asDense()` and `asPacked()` convert freely between the three representations; `asDense()` also promotes a plain array to a `denseImage`, without touching data already in one.


``` r
sparse <- asSparse(image > 2)
sparse
## Sparse image: 4 x 5 x 6 x 10 (logical)
##   Spatial dimensions : 4 x 5 x 6
##   Voxel size         : 1 x 1 x 1 (unit unknown)
##   Values per location: 10
##   Locations stored   : 29 of 120 (75.8% sparse)
sparseness(sparse)
## [1] 0.7583333
```

A sparse image's mask can be recovered as a logical array with `mask()`, and its stored values as a matrix with one column per present location with `maskedMatrix()`. Because a `sparseImage` stores its values already packed that way, this does not involve a copy.


``` r
mask(sparse)[1:2, 1:2, 1]
##       [,1] [,2]
## [1,]  TRUE TRUE
## [2,] FALSE TRUE
dim(maskedMatrix(sparse))
## [1] 10 29
```

`asPacked()` narrows storage to a narrow integer or `float32` type, choosing a slope and intercept automatically so the type's whole range is used.


``` r
packed <- asPacked(image, type = "int16")
packed
## Packed image: 4 x 5 x 6 x 10 (int16)
##   Spatial dimensions : 4 x 5 x 6
##   Voxel size         : 2 x 2 x 2.5 (unit unknown)
##   Values per location: 10
##   Scaling            : value = stored * 0.000107171 + 0.298598
##   Storage            : 2,400 bytes, against 9,600 as double
```

`asDense()` unpacks whichever compact representation it is given, so code that just wants ordinary values need not ask which one it has:


``` r
unpacked <- asDense(packed)
range(as.array(unpacked) - as.array(image)) # narrowing loses some precision, but very little
## [1] -5.334240e-05  5.350487e-05
```

`storageType()` reports `"double"`, `"integer"` and so on for anything that is not a packed image, and the packed type otherwise, so it can be used without first checking which kind of image is in hand. `isDenseImage()`, `isSparseImage()` and `isPackedImage()` (and the umbrella `isImage()`) are the best way to test for image types.


``` r
storageType(image)
## [1] "double"
storageType(packed)
## [1] "int16"
isPackedImage(packed)
## [1] TRUE
```

A `template` argument, accepted by all three constructors and by `asDense()`/`asSparse()`/`asPacked()`, carries geometry over from an existing image when only the data is changing.


``` r
doubled <- denseImage(as.array(image) * 2, template = image)
voxelSize(doubled)
## [1] 2.0 2.0 2.5
```

## Applying and reducing functions

`imapply()` is the memory-efficient analogue of `apply()` from the `base` package: it gives the same answer for the same `margin`, but gathers each sub-array directly through the image's strides rather than permuting the whole array into a fresh copy first.


``` r
totals <- imapply(image, 4, sum)
totals
##  [1]  13.1653065  -7.9258345  11.3782881  -6.8365127  -2.9111387 -18.8599395
##  [7]  -0.4747438   7.4987982 -12.3572488 -10.5496616
```

Three further verbs differ only in how much of the *space* is handed to the function at a time — a single location, a line, or a slice — with the values held at each location always travelling along with it. For the image above with dimensions 4 x 5 x 6 x 10 (the first three spatial), `voxelApply()` calls the function with 10 elements each time, `lineApply(... axis = 1)` passes a 4 x 10 matrix and `sliceApply(... axis = 3)` passes a 4 x 5 x 10 array.


``` r
means <- voxelApply(image, mean)
means # a single value per location comes back as an image, with the same geometry
## Dense image: 4 x 5 x 6 (double)
##   Spatial dimensions : 4 x 5 x 6
##   Voxel size         : 2 x 2 x 2.5 (unit unknown)

lineMeans <- lineApply(image, rowMeans, axis = 1)
dim(lineMeans)
## [1] 4 5 6

sliceMeans <- sliceApply(image, mean, axis = 3)
sliceMeans
## [1] -0.028966835 -0.074629673  0.013742921 -0.078449109 -0.006249744
## [6]  0.035189006
```

`voxelApply()` accepts a `mask` – a logical array, a sparse image, or any numeric array where non-zero elements are selected — and visits only those locations, leaving `fill` (`0` by default) everywhere else. The loop runs over a packed matrix of just the selected columns.


``` r
masked <- voxelApply(image, mean, mask = sparse)
masked[1:2, 1:2, 1]
##            [,1]        [,2]
## [1,] -0.1243711  0.05416769
## [2,]  0.0000000 -0.17414455
```

For a long-running call, `progress = TRUE` draws a text bar showing percentage complete and throughput in voxels per second. Alternatively, passing a function with arguments `done` and `total` allows you to report however you like.


``` r
voxelApply(image, slow_function, progress = TRUE)
##   |======================            | 62%  1.4M voxels/s
```

### Reductions

A handful of summary functions — `min`, `max`, `range`, `which.min`, `which.max`, plus `sum`, `mean`, `prod`, `var`, `sd`, `any`, `all` and `countNA` — can be computed without calling back into R at all, using `imreduce()`:


``` r
imreduce(image, 4, "range")
##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
## [1,] -2.214700 -2.888921 -2.403096 -3.008049 -2.528501 -2.939774 -2.596111
## [2,]  2.401618  2.497662  2.649167  2.165369  3.810277  2.675741  2.349493
##           [,8]      [,9]     [,10]
## [1,] -2.402231 -2.996949 -3.213189
## [2,]  3.055742  2.401222  2.185438
```

`imapply()` recognises calls to the base equivalents of `min`, `max`, `range`, `which.min` and `which.max`, and routes to `imreduce()` automatically, since these are the only ones guaranteed to give the same answer either way. Arithmetic reductions like `sum()` can differ from `imapply(x, margin, sum)` very slightly, because R accumulates `sum()` and `mean()` using the `long double` type internally.


``` r
identical(imapply(image, 4, max), imreduce(image, 4, "max"))
## [1] TRUE
```

## Parallelism

Work is parallelised in one of two ways, according to what is being run. Compiled reductions run concurrently in-process, using Grand Central Dispatch or OpenMP where available. The `parallelBackend()` function identifies which is in use in the current build.

An R function passed to `imapply()` and its variants cannot run on a worker thread — R's interpreter is not reentrant — so that work is instead parallelised by forking the R session using the `parallel` package, which is unavailable on Windows (`canFork()` reports whether it is available). Pass a `threads` argument to any
apply or reduce call, or set it globally via the `imply.threads` option:


``` r
options(imply.threads = 4)
```

## C++ API

The apply/reduce engine, the accessors that let one kernel body read dense, sparse and narrowly-stored images alike, and the voxel-to-world geometry machinery behind `worldTransform()` and friends, are all available to other packages' compiled code. Everything under `inst/include/` is header-only, so using it costs nothing to link against. The C++ headers include comments explaining their purpose and usage.

```
LinkingTo: imply
```

in a package's `DESCRIPTION` is enough; there is no library to link, and only the one header to include.

```c++
#include "imply.h"
```

The pieces most likely to be wanted from outside are as follows.

- **`Raster<D>`** — an *n*-dimensional index space with general strides, split at a runtime `spatial` index into leading spatial dimensions and trailing value dimensions. `FixedRaster<D>` keeps extents and indices on the stack when the dimensionality is known while compiling; `DynamicRaster` covers other cases.
- **`OffsetWalker`** — odometer traversal of an arbitrary subset of an image's dimensions, yielding memory offsets one addition at a time with no allocation. This is what lets a sub-array be gathered without permuting the whole image first.
- **`ImageSpace`** and **`Affine`** — voxel-to-world geometry, with the point conversions and rounding strategies that go with it, and no dependency on any file format.
- **`parallelFor`** — divides a range into chunks and runs them concurrently using libdispatch, OpenMP or a plain loop, whichever was compiled in. Nothing passed to it may touch the R API, allocate R objects, or draw from R's RNG, none of which are thread-safe.
- **`dispatchType`**, **`dispatchNarrowType`** and **`dispatchDims`** — the R-boundary switches from a runtime `SEXP`'s storage mode (or dimensionality) to a compile-time type tag, meant to be called exactly once per entry point so that everything below it is fully typed.
- **`DenseAccessor`**, **`SparseAccessor`** and **`NarrowAccessor`**, plus **`LocationMask`** — an `operator[]` over a plain pointer, a masked and packed array, or narrow storage, respectively, sharing one interface so that a kernel written against it serves all three without modification.
- **`Sink`**, **`VectorSink`**, **`TypedSink`** and **`ListSink`** — preallocated result accumulation for a per-call loop, falling back to a list only once a result turns out not to fit a plain vector.

A kernel that sums a numeric vector, written once against the shared accessor interface, therefore works whether the underlying storage is a plain R vector or something more compact. For example

```c++
// [[Rcpp::depends(imply)]]
#include "imply.h"
#include <Rcpp.h>

using namespace imply;

// One kernel body serves dense, sparse and narrow storage alike;
// only the accessor differs
template <typename Accessor>
double sumOver (const Accessor &values, const Extent n)
{
    double total = 0;
    for (Extent i = 0; i < n; i++)
        total += static_cast<double>(values[i]);
    return total;
}

// [[Rcpp::export]]
double denseSum (Rcpp::NumericVector x)
{
    DenseAccessor<double> accessor(REAL(x));
    return sumOver(accessor, Rf_xlength(x));
}
```

Here the templated `sumOver()` kernel function can handle any of the image representations, but only has to deal with indexing over an array of doubles. The `denseSum()` function shows an example of it being applied to a dense image.
