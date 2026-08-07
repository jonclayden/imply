#ifndef _IMPLY_PARALLEL_H_
#define _IMPLY_PARALLEL_H_

#include <algorithm>
#include <cstddef>
#include <string>

#if defined(HAVE_LIBDISPATCH)
#include <dispatch/dispatch.h>

// The system headers libdispatch pulls in define FALSE and TRUE as macros,
// which shadow the enumerators of R's Rboolean and make calls such as
// R_useDynamicSymbols(dll, FALSE) fail to resolve. R's own Boolean.h undefines
// them for the same reason, but only when it happens to be included later --
// and Rcpp::compileAttributes() puts a package's own header first
#undef FALSE
#undef TRUE

#elif defined(_OPENMP)
#include <omp.h>
#endif

namespace imply {

// Parallelism over three backends: libdispatch where available (which on
// macOS is the practical choice, since the system compiler rejects -fopenmp),
// OpenMP otherwise, and a plain loop when neither is present.
//
// Two departures from the equivalent shim in mmand:
//
//  * The work is expressed as a template taking a range, not as macros around
//    a loop body. mmand has to use macros because dispatch_apply() takes a
//    block, which is a language extension rather than a C++ construct;
//    dispatch_apply_f() takes an ordinary function pointer and a context, so
//    a capture-less lambda serves as the trampoline and no -fblocks support
//    or Blocks runtime is needed. That also removes the PARALLEL_LOOP_CONTINUE
//    leak, where `continue` had to become `return` depending on the backend.
//
//  * Work is divided into a fixed number of chunks, and the backend iterates
//    over chunks rather than over raw items. This is what makes a requested
//    thread count mean something under libdispatch, which offers no width
//    control of its own: with only n chunks to run, no more than n can be in
//    flight. In mmand the equivalent option is silently ignored whenever the
//    libdispatch path is taken.
//
// NOTHING passed to parallelFor may touch the R API, allocate R objects, or
// draw from R's global RNG, none of which are thread-safe.

inline const char * parallelBackend ()
{
#if defined(HAVE_LIBDISPATCH)
    return "libdispatch";
#elif defined(_OPENMP)
    return "openmp";
#else
    return "none";
#endif
}

inline bool parallelAvailable ()
{
#if defined(HAVE_LIBDISPATCH) || defined(_OPENMP)
    return true;
#else
    return false;
#endif
}

// The number of chunks to divide `items` into. A non-positive request means
// "as many as the backend sees fit", which is capped so that each chunk still
// holds at least `minimumPerChunk` items and the overhead stays worthwhile
inline std::size_t chunkCount (const std::size_t items, const int threads,
                               const std::size_t minimumPerChunk = 1)
{
    if (items == 0 || !parallelAvailable())
        return (items == 0 ? 0 : 1);

    std::size_t requested;
    if (threads > 0)
        requested = static_cast<std::size_t>(threads);
    else
    {
#if defined(_OPENMP)
        requested = static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
#else
        // libdispatch manages its own width, so ask for a reasonable number of
        // chunks and let it decide how many to run at once
        requested = 8;
#endif
    }

    const std::size_t affordable = std::max<std::size_t>(1, items / std::max<std::size_t>(1, minimumPerChunk));
    return std::max<std::size_t>(1, std::min(requested, std::min(items, affordable)));
}

namespace internal {

template <typename Functor>
struct ChunkContext
{
    Functor *fn;
    std::size_t items, chunks;

    void run (const std::size_t index) const
    {
        const std::size_t size = (items + chunks - 1) / chunks;
        const std::size_t begin = index * size;
        if (begin < items)
            (*fn)(begin, std::min(begin + size, items));
    }
};

} // namespace internal

// Split `items` into chunks and run `fn(begin, end)` on each, possibly
// concurrently. The ranges are disjoint and cover the whole extent, so a
// kernel writing only within its own range needs no synchronisation
template <typename Functor>
inline void parallelFor (const std::size_t items, const int threads, Functor fn)
{
    if (items == 0)
        return;

    const std::size_t chunks = chunkCount(items, threads);
    if (chunks <= 1)
    {
        fn(std::size_t(0), items);
        return;
    }

    internal::ChunkContext<Functor> context { &fn, items, chunks };

#if defined(HAVE_LIBDISPATCH)
    dispatch_apply_f(chunks, DISPATCH_APPLY_AUTO, &context, [](void *raw, std::size_t index) {
        static_cast<internal::ChunkContext<Functor> *>(raw)->run(index);
    });
#elif defined(_OPENMP)
    const long count = static_cast<long>(chunks);
    #pragma omp parallel for num_threads(count) schedule(static)
    for (long index=0; index<count; index++)
        context.run(static_cast<std::size_t>(index));
#else
    for (std::size_t index=0; index<chunks; index++)
        context.run(index);
#endif
}

} // namespace imply

#endif
