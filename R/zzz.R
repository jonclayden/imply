#' @keywords internal
#'
#' @useDynLib imply, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd var
#' @importFrom utils flush.console
"_PACKAGE"

## S7 methods defined at the top level of a package are recorded but not
## activated until the namespace is loaded, so they have to be registered here.
## Without this, print() and as.array() silently fall back to the array methods
.onLoad <- function (libname, pkgname)
{
    S7::S7_on_load()
    registerSparseMethods()
    registerPackedMethods()
    registerReductions()
}

.onUnload <- function (libpath)
{
    S7::S7_on_unload()
}

## method<- leaves a placeholder in the namespace for each generic owned by
## another package, such as print; this removes them once every top-level
## registration has run, which is why it has to stay in the last file collated
S7::S7_on_build()
