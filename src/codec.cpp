#include <Rcpp.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "imply/Raster.h"
#include "imply/Storage.h"
#include "imply/Dispatch.h"
#include "imply/Blocks.h"
#include "imply/View.h"

using namespace imply;

// Turning a byte stream into an image, and back. Nothing here knows about any
// file format: a format adapter parses its own header into a storage
// descriptor (type, byte order, scaling, layout) and hands over a connection
// positioned at the data. The bytes arrive through R connections, so that
// compression is R's business rather than ours, and are reinterpreted here.
//
// Colour types never reach this file. They are read as uint8 with an extra,
// fastest-varying storage axis, which the R side arranges through the layout
namespace {

enum class Code { int8, uint8, int16, uint16, int32, uint32, int64, float32, float64, complex64, complex128, bit };

Code codeFromName (const std::string &name)
{
    if (name == "int8")       return Code::int8;
    if (name == "uint8")      return Code::uint8;
    if (name == "int16")      return Code::int16;
    if (name == "uint16")     return Code::uint16;
    if (name == "int32")      return Code::int32;
    if (name == "uint32")     return Code::uint32;
    if (name == "int64")      return Code::int64;
    if (name == "float32")    return Code::float32;
    if (name == "float64")    return Code::float64;
    if (name == "complex64")  return Code::complex64;
    if (name == "complex128") return Code::complex128;
    if (name == "bit")        return Code::bit;
    Rcpp::stop("Unknown stored data type \"%s\"", name);
}

// 2^53: beyond this, not every integer has an exact double
const double exactLimit = 9007199254740992.0;

// Bytes are read unaligned, and reversed when the stream's byte order is not
// ours. A complex value is two components, each reversed in place
template <typename T>
inline T load (const Rbyte *p, const bool swap)
{
    T value;
    if (!swap)
        std::memcpy(&value, p, sizeof(T));
    else
    {
        Rbyte reversed[sizeof(T)];
        for (std::size_t i=0; i<sizeof(T); i++)
            reversed[i] = p[sizeof(T) - 1 - i];
        std::memcpy(&value, reversed, sizeof(T));
    }
    return value;
}

template <typename T>
inline void store (Rbyte *p, const T value, const bool swap)
{
    std::memcpy(p, &value, sizeof(T));
    if (swap)
        std::reverse(p, p + sizeof(T));
}

// Decodes one stored element into the R type the image will hold. Integer
// types stay integer unless they are scaled; wider integers and floats become
// double; complex stays complex. Scaling is not applied to complex data, as in
// NIfTI
template <typename Stored, typename OutTag, bool Complex = false>
struct Decoder
{
    typedef OutTag Tag;
    typedef typename OutTag::Type Out;
    static constexpr std::size_t size = (Complex ? 2 : 1) * sizeof(Stored);

    bool swap, scaled;
    double slope, intercept;

    Decoder (const bool swap, const double slope, const double intercept)
        : swap(swap), scaled(slope != 1.0 || intercept != 0.0), slope(slope), intercept(intercept) {}

    Out operator() (const Rbyte *p) const
    {
        if constexpr (Complex)
        {
            Rcomplex z;
            z.r = static_cast<double>(load<Stored>(p, swap));
            z.i = static_cast<double>(load<Stored>(p + sizeof(Stored), swap));
            return z;
        }
        else if constexpr (OutTag::kind == StorageType::integer)
        {
            const Stored value = load<Stored>(p, swap);
            if constexpr (std::is_same<Stored, std::int32_t>::value)
            {
                // R reserves the smallest int for NA
                if (value == INT_MIN)
                    Rcpp::stop("The stored value %d cannot be represented as an R integer; read the image as packed instead", INT_MIN);
            }
            return static_cast<int>(value);
        }
        else
        {
            const Stored stored = load<Stored>(p, swap);
            if constexpr (std::is_same<Stored, float>::value)
            {
                // As for narrow storage: a float NaN stands for a missing value
                if (stored != stored)
                    return NA_REAL;
            }
            if constexpr (std::is_same<Stored, std::int64_t>::value)
            {
                if (std::fabs(static_cast<double>(stored)) > exactLimit)
                    Rcpp::stop("A stored 64-bit integer exceeds 2^53, and cannot be represented exactly as a double");
            }
            const double value = static_cast<double>(stored);
            return scaled ? value * slope + intercept : value;
        }
    }
};

// Run fn with the decoder for a stored type, chosen once here so that the
// loops below are fully typed
template <class Functor>
SEXP withDecoder (const Code code, const bool swap, const double slope, const double intercept, Functor &&fn)
{
    const bool scaled = (slope != 1.0 || intercept != 0.0);
    switch (code)
    {
        case Code::int8:
        if (scaled) return fn(Decoder<std::int8_t, RealTag>(swap, slope, intercept));
        else        return fn(Decoder<std::int8_t, IntegerTag>(swap, slope, intercept));

        case Code::uint8:
        if (scaled) return fn(Decoder<std::uint8_t, RealTag>(swap, slope, intercept));
        else        return fn(Decoder<std::uint8_t, IntegerTag>(swap, slope, intercept));

        case Code::int16:
        if (scaled) return fn(Decoder<std::int16_t, RealTag>(swap, slope, intercept));
        else        return fn(Decoder<std::int16_t, IntegerTag>(swap, slope, intercept));

        case Code::uint16:
        if (scaled) return fn(Decoder<std::uint16_t, RealTag>(swap, slope, intercept));
        else        return fn(Decoder<std::uint16_t, IntegerTag>(swap, slope, intercept));

        case Code::int32:
        if (scaled) return fn(Decoder<std::int32_t, RealTag>(swap, slope, intercept));
        else        return fn(Decoder<std::int32_t, IntegerTag>(swap, slope, intercept));

        case Code::uint32:     return fn(Decoder<std::uint32_t, RealTag>(swap, slope, intercept));
        case Code::int64:      return fn(Decoder<std::int64_t, RealTag>(swap, slope, intercept));
        case Code::float32:    return fn(Decoder<float, RealTag>(swap, slope, intercept));
        case Code::float64:    return fn(Decoder<double, RealTag>(swap, slope, intercept));
        case Code::complex64:  return fn(Decoder<float, ComplexTag, true>(swap, slope, intercept));
        case Code::complex128: return fn(Decoder<double, ComplexTag, true>(swap, slope, intercept));
        case Code::bit:        break;
    }
    Rcpp::stop("Bit data are decoded separately");
}

// Pulls bytes through R. The callbacks are R closures around the connection;
// Rcpp::Function evaluates through Rcpp_fast_eval, so an error or interrupt in
// R unwinds this stack cleanly
class Source
{
protected:
    Rcpp::Function read_, skip_;

public:
    // Bytes are read in chunks of at most this many, which bounds the
    // transient memory for a large image and keeps each read within the
    // length of a standard R vector
    static constexpr Extent chunkBytes = 1u << 26;

    Source (Rcpp::Function read, Rcpp::Function skip) : read_(read), skip_(skip) {}

    Rcpp::RawVector read (const Extent n)
    {
        Rcpp::RawVector bytes = read_(static_cast<double>(n));
        if (static_cast<Extent>(bytes.size()) != n)
            Rcpp::stop("Unexpected end of data: %.0f bytes were wanted but only %.0f were available",
                       static_cast<double>(n), static_cast<double>(bytes.size()));
        return bytes;
    }

    void skip (const Extent n)
    {
        if (n > 0)
            skip_(static_cast<double>(n));
    }

    // Call fn on each of `count` consecutive elements of `size` bytes, in
    // stream order, reading as it goes
    template <class Functor>
    void forEach (const Extent count, const std::size_t size, Functor &&fn)
    {
        const Extent perChunk = std::max<Extent>(1, chunkBytes / size);
        for (Extent done=0; done<count; )
        {
            const Extent n = std::min(perChunk, count - done);
            const Rcpp::RawVector bytes = read(n * size);
            const Rbyte *p = bytes.begin();
            for (Extent j=0; j<n; j++)
                fn(p + j * size);
            done += n;
        }
    }

    void readInto (Rbyte *out, const Extent n)
    {
        for (Extent done=0; done<n; )
        {
            const Extent m = std::min(chunkBytes, n - done);
            const Rcpp::RawVector bytes = read(m);
            std::copy(bytes.begin(), bytes.end(), out + done);
            done += m;
        }
    }
};

// Blocks requested, each with the positions in the result it fills. Reading
// them in stream order means a compressed stream is only ever read forwards,
// and a block asked for twice is read once
struct BlockPlan
{
    std::vector<std::pair<Extent, std::vector<Extent>>> reads;
    Extent slots;
};

BlockPlan planBlocks (const Rcpp::NumericVector &blocks)
{
    std::vector<std::pair<Extent, Extent>> wanted;
    for (R_xlen_t k=0; k<blocks.size(); k++)
        wanted.push_back(std::make_pair(static_cast<Extent>(blocks[k]), static_cast<Extent>(k)));
    std::sort(wanted.begin(), wanted.end());

    BlockPlan plan;
    plan.slots = static_cast<Extent>(blocks.size());
    for (std::size_t i=0; i<wanted.size(); i++)
    {
        if (plan.reads.empty() || plan.reads.back().first != wanted[i].first)
            plan.reads.push_back(std::make_pair(wanted[i].first, std::vector<Extent>()));
        plan.reads.back().second.push_back(wanted[i].second);
    }
    return plan;
}

// The whole stream, in any layout, into a dense array in view order. The
// stream is visited in its own order, and the inverse walker says where each
// element belongs in the view
template <class D>
SEXP decodeDense (Source &source, const D &decode, const ViewMap &view)
{
    typedef typename D::Tag Tag;
    Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(view.size()));
    typename Tag::Type *out = result.begin();

    OffsetWalker walker = view.inverseWalker();
    walker.reset();
    source.forEach(view.size(), D::size, [&](const Rbyte *p) {
        out[view.inverseBase + walker.offset()] = decode(p);
        walker.next();
    });
    return result;
}

// Selected volumes of block-ordered data, into a dense array whose last
// dimension runs over the volumes asked for. Within a block only the spatial
// axes can be reordered, which `spatial` describes
template <class D>
SEXP decodeDenseBlocks (Source &source, const D &decode, const ViewMap &spatial, const BlockPlan &plan)
{
    typedef typename D::Tag Tag;
    const Extent blockElements = spatial.size();
    const Extent blockBytes = blockElements * D::size;
    Rcpp::Vector<Tag::sexpType> result(static_cast<R_xlen_t>(blockElements * plan.slots));
    typename Tag::Type *out = result.begin();

    Extent position = 0;
    for (std::size_t r=0; r<plan.reads.size(); r++)
    {
        const Extent block = plan.reads[r].first;
        const std::vector<Extent> &slots = plan.reads[r].second;
        source.skip((block - position) * blockBytes);

        typename Tag::Type *first = out + slots[0] * blockElements;
        OffsetWalker walker = spatial.inverseWalker();
        walker.reset();
        source.forEach(blockElements, D::size, [&](const Rbyte *p) {
            first[spatial.inverseBase + walker.offset()] = decode(p);
            walker.next();
        });
        for (std::size_t k=1; k<slots.size(); k++)
            std::copy(first, first + blockElements, out + slots[k] * blockElements);

        position = block + 1;
    }
    return result;
}

std::size_t maskBytes (const Extent locations) { return ((locations + 63) / 64) * 8; }

inline void setBit (Rbyte *bytes, const Extent i) { bytes[i >> 3] |= static_cast<Rbyte>(1u << (i & 7)); }

inline bool getBit (const Rbyte *bytes, const Extent i) { return (bytes[i >> 3] >> (i & 7)) & 1u; }

// Selected volumes of block-ordered data, into the parts of a sparse image:
// a mask over locations in stored order, and a matrix of values with one
// column per stored location. With a mask given, exactly the locations it
// selects are kept, and the result is built in one pass. Without one, a
// location is kept if any value there is not zero, which needs the non-zero
// values set aside until every block has been seen
template <class D>
Rcpp::List decodeSparseBlocks (Source &source, const D &decode, const ViewMap &spatial,
                               const BlockPlan &plan, SEXP viewMask)
{
    typedef typename D::Tag Tag;
    typedef typename Tag::Type Value;
    const Extent locations = spatial.size();
    const Extent blockBytes = locations * D::size;
    const Extent slots = plan.slots;

    Rcpp::RawVector mask(static_cast<R_xlen_t>(maskBytes(locations)));
    std::fill(mask.begin(), mask.end(), Rbyte(0));
    Rbyte * const bits = mask.begin();
    std::vector<Extent> rank;

    // Each location's place among those stored
    auto rankAll = [&]() {
        rank.assign(locations, 0);
        Extent count = 0;
        for (Extent s=0; s<locations; s++)
        {
            if (getBit(bits, s))
                rank[s] = count++;
        }
        return count;
    };

    if (!Rf_isNull(viewMask))
    {
        // The mask is in view order; walking the storage order with the
        // inverse map says which stored location each view location is
        const Rcpp::LogicalVector selected(viewMask);
        OffsetWalker walker = spatial.inverseWalker();
        walker.reset();
        for (Extent s=0; s<locations; s++)
        {
            if (selected[spatial.inverseBase + walker.offset()])
                setBit(bits, s);
            walker.next();
        }

        const Extent stored = rankAll();
        Rcpp::Vector<Tag::sexpType> values(static_cast<R_xlen_t>(slots * stored));
        Value *out = values.begin();

        Extent position = 0;
        for (std::size_t r=0; r<plan.reads.size(); r++)
        {
            const Extent block = plan.reads[r].first;
            const std::vector<Extent> &targets = plan.reads[r].second;
            source.skip((block - position) * blockBytes);

            Extent s = 0;
            source.forEach(locations, D::size, [&](const Rbyte *p) {
                if (getBit(bits, s))
                {
                    const Value value = decode(p);
                    for (std::size_t k=0; k<targets.size(); k++)
                        out[rank[s] * slots + targets[k]] = value;
                }
                s++;
            });
            position = block + 1;
        }

        values.attr("dim") = Rcpp::Dimension(static_cast<int>(slots), static_cast<int>(stored));
        return Rcpp::List::create(Rcpp::Named("mask") = mask, Rcpp::Named("values") = values);
    }

    std::vector<Extent> where, slot;
    std::vector<Value> kept;

    Extent position = 0;
    for (std::size_t r=0; r<plan.reads.size(); r++)
    {
        const Extent block = plan.reads[r].first;
        const std::vector<Extent> &targets = plan.reads[r].second;
        source.skip((block - position) * blockBytes);

        Extent s = 0;
        source.forEach(locations, D::size, [&](const Rbyte *p) {
            const Value value = decode(p);
            if (!Tag::isZero(value))
            {
                setBit(bits, s);
                for (std::size_t k=0; k<targets.size(); k++)
                {
                    where.push_back(s);
                    slot.push_back(targets[k]);
                    kept.push_back(value);
                }
            }
            s++;
        });
        position = block + 1;
    }

    const Extent stored = rankAll();
    Rcpp::Vector<Tag::sexpType> values(static_cast<R_xlen_t>(slots * stored));
    std::fill(values.begin(), values.end(), Value());
    Value *out = values.begin();
    for (std::size_t i=0; i<kept.size(); i++)
        out[rank[where[i]] * slots + slot[i]] = kept[i];

    values.attr("dim") = Rcpp::Dimension(static_cast<int>(slots), static_cast<int>(stored));
    return Rcpp::List::create(Rcpp::Named("mask") = mask, Rcpp::Named("values") = values);
}

// Bit data are small enough to read whole. Element e of the stream is bit
// e % 8 of byte e / 8, counting from the least or most significant end
inline bool bitAt (const Rbyte *bytes, const Extent e, const bool msbFirst)
{
    const unsigned shift = static_cast<unsigned>(e & 7);
    return (bytes[e >> 3] >> (msbFirst ? 7 - shift : shift)) & 1u;
}

SEXP decodeBits (Source &source, const ViewMap &view, SEXP blocks, const ViewMap &spatial, const Extent total,
                 const bool msbFirst)
{
    const Rcpp::RawVector bytes = source.read((total + 7) / 8);
    const Rbyte *data = bytes.begin();

    if (Rf_isNull(blocks))
    {
        Rcpp::LogicalVector result(static_cast<R_xlen_t>(view.size()));
        OffsetWalker walker = view.inverseWalker();
        walker.reset();
        for (Extent e=0; e<view.size(); e++)
        {
            result[view.inverseBase + walker.offset()] = bitAt(data, e, msbFirst);
            walker.next();
        }
        return result;
    }

    const Rcpp::NumericVector wanted(blocks);
    const Extent blockElements = spatial.size();
    Rcpp::LogicalVector result(static_cast<R_xlen_t>(blockElements * wanted.size()));
    for (R_xlen_t k=0; k<wanted.size(); k++)
    {
        // A block need not start on a byte boundary, so the bit offset is
        // carried through rather than rounded
        const Extent start = static_cast<Extent>(wanted[k]) * blockElements;
        OffsetWalker walker = spatial.inverseWalker();
        walker.reset();
        for (Extent e=0; e<blockElements; e++)
        {
            result[k * blockElements + spatial.inverseBase + walker.offset()] = bitAt(data, start + e, msbFirst);
            walker.next();
        }
    }
    return result;
}

void swapInPlace (Rbyte *bytes, const Extent n, const std::size_t unit)
{
    if (unit < 2)
        return;
    for (Extent i=0; i<n; i+=unit)
        std::reverse(bytes + i, bytes + i + unit);
}

std::vector<Extent> extentsOf (const Rcpp::IntegerVector &dim)
{
    for (R_xlen_t i=0; i<dim.size(); i++)
    {
        if (dim[i] == NA_INTEGER || dim[i] < 0)
            Rcpp::stop("Dimensions must not be missing or negative");
    }
    return std::vector<Extent>(dim.begin(), dim.end());
}

} // anonymous namespace

// Decode a stream into dense values, the bytes of a packed image, or the
// parts of a sparse image. `dim` and `layout` describe the whole stream; when
// `blocks` is given, only those volumes (zero-based, in the order wanted) are
// read, and the stream must be stored a volume at a time. The connection is
// already positioned at the first byte of data
// [[Rcpp::export]]
SEXP decodeImage (Rcpp::Function read, Rcpp::Function skip, std::string type, bool swap,
                  double slope, double intercept, Rcpp::IntegerVector dim, SEXP layout,
                  int spatial, std::string as, SEXP blocks = R_NilValue, SEXP mask = R_NilValue,
                  bool msbFirst = false)
{
    Source source(read, skip);
    const Code code = codeFromName(type);
    const std::vector<Extent> dims = extentsOf(dim);
    const std::vector<int> order = layoutFrom(layout);
    const ViewMap view(dims, order);

    if (spatial < 0 || spatial > static_cast<int>(dims.size()))
        Rcpp::stop("The number of spatial dimensions is out of range");

    // The spatial part of a block-ordered layout, which is a layout in its
    // own right when blocks are being read. Other layouts interleave the
    // spatial axes with the rest, and are only ever read whole
    bool blockOrdered = true;
    std::vector<int> spatialOrder;
    if (!order.empty())
    {
        for (int i=0; i<static_cast<int>(order.size()); i++)
        {
            const int rank = std::abs(order[i]);
            if ((i < spatial && rank > spatial) || (i >= spatial && order[i] != i + 1))
                blockOrdered = false;
        }
        if (blockOrdered)
            spatialOrder.assign(order.begin(), order.begin() + spatial);
    }
    const std::vector<Extent> spatialDims(dims.begin(), dims.begin() + spatial);
    const ViewMap spatialView(spatialDims, spatialOrder);

    if (!blockOrdered && (!Rf_isNull(blocks) || as == "sparse"))
        Rcpp::stop("Volumes can only be read separately from data stored a volume at a time");

    Rcpp::NumericVector blockList;
    if (!Rf_isNull(blocks))
        blockList = Rcpp::NumericVector(blocks);
    else if (as == "sparse")
    {
        // A sparse read is always block by block, so reading everything is
        // reading every block in turn
        const Extent n = view.size() / std::max<Extent>(1, spatialView.size());
        blockList = Rcpp::NumericVector(static_cast<R_xlen_t>(n));
        for (Extent k=0; k<n; k++)
            blockList[k] = static_cast<double>(k);
    }

    if (code == Code::bit)
    {
        if (as != "dense")
            Rcpp::stop("Bit data can only be read densely");
        return decodeBits(source, view, blocks, spatialView, view.size(), msbFirst);
    }

    if (as == "packed")
    {
        std::size_t size = 0;
        switch (code)
        {
            case Code::int8: case Code::uint8: size = 1; break;
            case Code::int16: case Code::uint16: size = 2; break;
            case Code::int32: case Code::float32: size = 4; break;
            default: Rcpp::stop("Data of type %s cannot be kept packed; read them densely instead", type);
        }

        if (Rf_isNull(blocks))
        {
            Rcpp::RawVector result(static_cast<R_xlen_t>(view.size() * size));
            source.readInto(result.begin(), view.size() * size);
            if (swap)
                swapInPlace(result.begin(), view.size() * size, size);
            return result;
        }

        const BlockPlan plan = planBlocks(blockList);
        const Extent blockBytes = spatialView.size() * size;
        Rcpp::RawVector result(static_cast<R_xlen_t>(blockBytes * plan.slots));
        Extent position = 0;
        for (std::size_t r=0; r<plan.reads.size(); r++)
        {
            const Extent block = plan.reads[r].first;
            const std::vector<Extent> &slots = plan.reads[r].second;
            source.skip((block - position) * blockBytes);
            Rbyte *first = result.begin() + slots[0] * blockBytes;
            source.readInto(first, blockBytes);
            if (swap)
                swapInPlace(first, blockBytes, size);
            for (std::size_t k=1; k<slots.size(); k++)
                std::copy(first, first + blockBytes, result.begin() + slots[k] * blockBytes);
            position = block + 1;
        }
        return result;
    }

    return withDecoder(code, swap, slope, intercept, [&](auto decode) -> SEXP {
        if (as == "sparse")
            return decodeSparseBlocks(source, decode, spatialView, planBlocks(blockList), mask);
        if (Rf_isNull(blocks))
            return decodeDense(source, decode, view);
        return decodeDenseBlocks(source, decode, spatialView, planBlocks(blockList));
    });
}

// Encode values given in view order into a stream laid out as described. Out
// of range and missing values are refused rather than clamped or dropped: a
// file silently holding something other than what was written is worse than
// an error
// [[Rcpp::export]]
Rcpp::RawVector encodeImage (Rcpp::RObject x, std::string type, bool swap, double slope, double intercept,
                             Rcpp::IntegerVector dim, SEXP layout = R_NilValue, bool msbFirst = false)
{
    const Code code = codeFromName(type);
    const ViewMap view(extentsOf(dim), layoutFrom(layout));
    const Extent n = view.size();
    if (static_cast<Extent>(Rf_xlength(x)) != n)
        Rcpp::stop("Dimensions imply %.0f values, but %.0f were given", double(n), double(Rf_xlength(x)));
    if (slope == 0.0 || ISNAN(slope) || ISNAN(intercept))
        Rcpp::stop("Scaling must be finite, with a non-zero slope");

    OffsetWalker walker = view.inverseWalker();

    if (code == Code::bit)
    {
        Rcpp::RawVector result(static_cast<R_xlen_t>((n + 7) / 8));
        std::fill(result.begin(), result.end(), Rbyte(0));
        const Rcpp::NumericVector values(Rcpp::as<Rcpp::NumericVector>(x));
        walker.reset();
        for (Extent e=0; e<n; e++)
        {
            const double value = values[view.inverseBase + walker.offset()];
            if (ISNAN(value))
                Rcpp::stop("Missing values cannot be stored as bits");
            if (value != 0.0)
            {
                const unsigned shift = static_cast<unsigned>(e & 7);
                result[e >> 3] |= static_cast<Rbyte>(1u << (msbFirst ? 7 - shift : shift));
            }
            walker.next();
        }
        return result;
    }

    if (code == Code::complex64 || code == Code::complex128)
    {
        const Rcpp::ComplexVector values(Rcpp::as<Rcpp::ComplexVector>(x));
        const std::size_t part = (code == Code::complex64 ? 4 : 8);
        Rcpp::RawVector result(static_cast<R_xlen_t>(n * 2 * part));
        Rbyte *out = result.begin();
        walker.reset();
        for (Extent e=0; e<n; e++)
        {
            const Rcomplex z = values[view.inverseBase + walker.offset()];
            if (part == 4)
            {
                store<float>(out + e * 8, static_cast<float>(z.r), swap);
                store<float>(out + e * 8 + 4, static_cast<float>(z.i), swap);
            }
            else
            {
                store<double>(out + e * 16, z.r, swap);
                store<double>(out + e * 16 + 8, z.i, swap);
            }
            walker.next();
        }
        return result;
    }

    if (TYPEOF(x) == CPLXSXP)
        Rcpp::stop("Complex values can only be stored as complex64 or complex128");
    const Rcpp::NumericVector values(Rcpp::as<Rcpp::NumericVector>(x));

    auto encode = [&](auto stored, const double low, const double high, const bool integral) -> Rcpp::RawVector {
        typedef decltype(stored) Stored;
        Rcpp::RawVector result(static_cast<R_xlen_t>(n * sizeof(Stored)));
        Rbyte *out = result.begin();
        walker.reset();
        for (Extent e=0; e<n; e++)
        {
            const double raw = values[view.inverseBase + walker.offset()];
            walker.next();

            if (ISNAN(raw))
            {
                if (integral)
                    Rcpp::stop("Missing values cannot be stored in type %s; use a floating-point type", type);
                // R's NA is a particular NaN; float32 keeps only that it is
                // not a number, while float64 keeps the value exactly
                store<Stored>(out + e * sizeof(Stored), static_cast<Stored>(raw), swap);
                continue;
            }

            double value = (raw - intercept) / slope;
            if (integral)
            {
                value = std::nearbyint(value);
                if (value < low || value > high)
                    Rcpp::stop("The value %g is outside the range of type %s under the scaling given", raw, type);
            }
            store<Stored>(out + e * sizeof(Stored), static_cast<Stored>(value), swap);
        }
        return result;
    };

    switch (code)
    {
        case Code::int8:    return encode(std::int8_t(), -128.0, 127.0, true);
        case Code::uint8:   return encode(std::uint8_t(), 0.0, 255.0, true);
        case Code::int16:   return encode(std::int16_t(), -32768.0, 32767.0, true);
        case Code::uint16:  return encode(std::uint16_t(), 0.0, 65535.0, true);
        case Code::int32:   return encode(std::int32_t(), -2147483648.0, 2147483647.0, true);
        case Code::uint32:  return encode(std::uint32_t(), 0.0, 4294967295.0, true);
        case Code::int64:   return encode(std::int64_t(), -exactLimit, exactLimit, true);
        case Code::float32: return encode(float(), R_NegInf, R_PosInf, false);
        case Code::float64: return encode(double(), R_NegInf, R_PosInf, false);
        default: break;
    }
    Rcpp::stop("Unhandled stored data type");
}

// Reverse the bytes of each element, for packed data going to a stream of the
// other byte order
// [[Rcpp::export]]
Rcpp::RawVector swapBytes (Rcpp::RawVector bytes, int size)
{
    Rcpp::RawVector result = Rcpp::clone(bytes);
    swapInPlace(result.begin(), static_cast<Extent>(result.size()), static_cast<std::size_t>(size));
    return result;
}
