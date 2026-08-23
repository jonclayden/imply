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
// Work is divided into a fixed number of chunks, and the backend iterates over
// chunks rather than over raw items. This is what makes a requested thread
// count mean something under libdispatch, which offers no width control of its
// own: with only n chunks to run, no more than n can be in flight.
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
// serial, not "as many as the backend feels like": a shared library has no
// way to know whether the process already owns a slice of a larger
// allocation (a forked worker, an HPC job, a server handling concurrent
// requests), so silently claiming every core would be the wrong default.
// Multi-core use is opt-in, via an explicit thread count
inline std::size_t chunkCount (const std::size_t items, const int threads,
                               const std::size_t minimumPerChunk = 1)
{
    if (items == 0 || !parallelAvailable())
        return (items == 0 ? 0 : 1);

    const std::size_t requested = (threads > 0) ? static_cast<std::size_t>(threads) : 1;
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
