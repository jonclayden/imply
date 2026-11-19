This file documents the significant user-visible changes in each release of the
`imply` R package.

# imply 0.2.0

## Breaking changes

- Image geometry is now a class of its own, `imageGeometry`, held by every
  image in a single `geometry` property. The `spatial`, `voxelSize`,
  `orientation`, `spaceUnit` and `timeUnit` properties are gone: use the
  accessor functions, such as `spatial()` and `voxelSize()`, or `geometry()`.
- The constructors' `template` argument is replaced by `geometry`, which takes
  an `imageGeometry` or an image to copy one from. Explicit geometry arguments
  still take precedence over it.
- `spaceUnit` and `timeUnit` are replaced by a single `unit`, for voxel size
  and world coordinates. There is no longer a time unit, or any other
  description of the value dimensions: what the values at a location mean is
  for client packages to model.
- `voxelApply()`, `lineApply()` and `sliceApply()` now return a dense image
  when the function gives one value per location, line or slice. Its geometry
  is that of the input with the traversed axes removed, each location sitting
  at the centre of what it summarises, so `centre()` is unchanged. A plain
  array in still gives a plain array out, and `imapply()` and `imreduce()` are
  unaffected.
- The package now requires S7 1.0.0 or later, which in turn requires R 4.2.0 or
  later. It no longer imports `methods`.

## New features

- New geometry functions: `geometry()` and `geometry<-`, `imageGeometry()`
  and `isImageGeometry()`, `centre()` (also spelled `center()`), `extent()`
  and `sameGeometry()`. Every geometry function accepts an image, a bare
  geometry or a plain array.
- Packed and sparse images can now be reoriented without copying any data.
  They carry a `layout` saying how their axes map onto storage, so that
  `reorient()` changes only the view; indexing, `dim()`, the geometry and
  every apply function work through the view. `worldAxes()` reports how an
  image's axes run relative to world space, and `storageLayout()` gives its
  layout. Dense images are reoriented by copying. No anatomical or other domain
  convention is assumed about the world axes.
- The new `storedValues()` function gives the values of a sparse image exactly
  as stored, without copying.
- New `storageDescriptor()`, `readImageData()` and `writeImageData()` functions
  read and write image data from and to R connections, for a range of storage
  types, byte orders, offsets, scalings and axis layouts, including bit data.
  They know nothing about file formats: a format package parses a header into
  a descriptor and hands over a connection. Sparse images are read a volume at
  a time where the data are stored volume by volume.
- Packed and sparse images now support index assignment (`[<-`) without
  materialising the data. Packed assignment keeps the existing scaling and
  refuses a value it cannot represent.
- Unary `+`, `-` and `!` now work on sparse and packed images. Negating a
  sparse image keeps it sparse.
- `voxelApply()` accepts a numeric array or any image as `mask`, as well as a
  logical array. `asMaskVector()`, which interprets these, is now exported.

## Other changes

- When choosing a scaling for, or writing to, int32 storage, the most
  negative integer is no longer used, since R reads it as `NA`.
- The C++ API (`IMPLY_API_VERSION` 2) adds `ViewMap`, for mapping an image's
  axes onto its storage, and `StrideIterator`, a random-access iterator over a
  strided run of values for use with standard algorithms.

# imply 0.1.0

- First public release.
