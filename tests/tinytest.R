## CRAN permits at most two cores during checks
options(imply.threads = 2L)

if (requireNamespace("tinytest", quietly=TRUE))
    tinytest::test_package("imply")
