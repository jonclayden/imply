## The public C++ header interface.

## Everything is header-only, so a package linking to imply needs no library.
## These checks are on the shipped layout rather than on behaviour
headers <- system.file("include", package = "imply")
expect_true(nzchar(headers))
expect_true(file.exists(file.path(headers, "imply.h")))

for (header in c("Raster.h", "Space.h", "Blocks.h", "Parallel.h", "Storage.h",
                 "Dispatch.h", "RImage.h", "Sparse.h", "Narrow.h", "Sink.h"))
    expect_true(file.exists(file.path(headers, "imply", header)), info = header)

## The umbrella pulls in every one of them
umbrella <- readLines(file.path(headers, "imply.h"))
for (header in c("Raster.h", "Space.h", "Blocks.h", "Parallel.h", "Sparse.h", "Narrow.h"))
    expect_true(any(grepl(paste0('include "imply/', header), umbrella, fixed = TRUE)), info = header)

## Space.h carries its own definitions, so there is nothing to link against
space <- readLines(file.path(headers, "imply", "Space.h"))
expect_true(any(grepl("^inline .*Affine::inverse", space)))

## Raster.h is free of any dependency on R, so it can be used on its own
raster <- readLines(file.path(headers, "imply", "Raster.h"))
expect_false(any(grepl("Rcpp\\.h|Rinternals\\.h|<R\\.h>", raster)))
