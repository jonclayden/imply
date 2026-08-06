#' @keywords internal
#'
#' @useDynLib imply, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods setMethod
"_PACKAGE"

## S7 methods defined at the top level of a package are recorded but not
## activated until the namespace is loaded, so they have to be registered here.
## Without this, print() and as.array() silently fall back to the array methods
.onLoad <- function (libname, pkgname)
{
    S7::methods_register()
    registerSparseMethods()
    registerPackedMethods()
}
