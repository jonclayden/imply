#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([imply], [0.1.0], [code@clayden.org])

AC_LANG([C++])

# Get compiler flags from R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_PROG_CXX

# Both backends are probed; which one is used is decided by the preprocessor
# in src/Parallel.h, where libdispatch takes precedence. On macOS the system
# compiler generally rejects -fopenmp, which is exactly why the libdispatch
# path exists
m4_include([tools/m4/ax_openmp.m4])
AX_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

m4_include([tools/m4/jc_dispatch.m4])
JC_DISPATCH
AC_SUBST(LIBDISPATCH_CPPFLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
