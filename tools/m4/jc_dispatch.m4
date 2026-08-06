# Detect support for libdispatch (aka Grand Central Dispatch).
#
# Adapted from the version in the mmand package, which additionally probes for
# the block syntax extension and the Blocks runtime. Those are only needed by
# dispatch_apply(), which takes a block; dispatch_apply_f() takes an ordinary
# function pointer and a context argument instead, so none of that machinery
# is required here.

AC_DEFUN([JC_DISPATCH], [

AC_SEARCH_LIBS([dispatch_apply_f], [dispatch])
AC_CHECK_HEADER([dispatch/dispatch.h])

AC_MSG_CHECKING([whether a simple program can be compiled against libdispatch])
AC_LINK_IFELSE([AC_LANG_SOURCE([[

#include <stdlib.h>
#include <dispatch/dispatch.h>

static void kernel (void *context, size_t iteration)
{
    int *values = (int *) context;
    values[iteration] = iteration;
}

int main ()
{
    int *values = (int *) calloc(10, sizeof(int));
    dispatch_apply_f(10, DISPATCH_APPLY_AUTO, values, &kernel);
    free(values);
    return 0;
}

]])], [
    LIBDISPATCH_CPPFLAGS="-DHAVE_LIBDISPATCH"
    AC_MSG_RESULT([yes])
], [AC_MSG_RESULT([no])])

])
