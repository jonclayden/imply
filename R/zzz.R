#' @keywords internal
#'
#' @useDynLib imply, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods setMethod
#' @importFrom stats sd var
#' @importFrom utils flush.console
"_PACKAGE"

## S7 methods defined at the top level of a package are recorded but not
## activated until the namespace is loaded, so they have to be registered here.
## Without this, print() and as.array() silently fall back to the array methods
.onLoad <- function (libname, pkgname)
{
    if ("S7_on_load" %in% getNamespaceExports("S7")) {
      getExportedValue("S7", "S7_on_load")()
    } else {
      S7::methods_register()
    }
    registerSparseMethods()
    registerPackedMethods()
    registerReductions()
}

.onUnload <- function(libpath) {
  if ("S7_on_unload" %in% getNamespaceExports("S7")) {
    getExportedValue("S7", "S7_on_unload")()
  }
}

if ("S7_on_build" %in% getNamespaceExports("S7")) {
  getExportedValue("S7", "S7_on_build")()
}
