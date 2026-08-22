#ifndef _IMPLY_SPACE_H_
#define _IMPLY_SPACE_H_

#include <Rcpp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace imply {

// A self-contained Mersenne twister, so that anything needing randomness owns
// its own stream rather than reaching for R's global RNG. R's generator is
// shared mutable state and cannot be touched from a worker thread; this can.
//
// Seeding it from R's RNG at the boundary keeps set.seed() in charge of
// reproducibility. To keep results independent of how work is divided between
// threads, a parallel caller should seed one generator per fixed-size block
// from the block's position, rather than one per thread
class RandomGenerator
{
protected:
    std::mt19937_64 engine;
    std::uniform_real_distribution<double> distribution;

public:
    explicit RandomGenerator (const std::uint64_t seed)
        : engine(seed), distribution(0.0, 1.0) {}

    // Uniform on [0,1), matching the convention of R's unif_rand()
    double uniform () { return distribution(engine); }

    void reseed (const std::uint64_t seed) { engine.seed(seed); }
};

// Location conventions: voxel-indexed, scaled for voxel dimensions only (as
// with a diagonal xform), or world coordinates fully respecting the xform
enum class PointType { voxel, scaled, world };

// Rounding strategies: none, standard for nearest-neighbour, or probabilistic
// for stochastic nearest neighbour (probabilities proportional to distance)
enum class RoundingType { none, conventional, probabilistic };

typedef std::array<double,3> Point;

// A 4x4 affine transform, stored row-major. Only the affine case is supported,
// meaning the final row is implicitly (0,0,0,1), which is what makes the
// inverse cheap and exact
class Affine
{
protected:
    std::array<double,16> values;

public:
    Affine () { values.fill(0.0); }

    static Affine identity ()
    {
        Affine result;
        for (int i=0; i<4; i++)
            result(i,i) = 1.0;
        return result;
    }

    // A diagonal transform built from voxel dimensions, used when an image
    // carries no explicit transform of its own
    static Affine scaling (const std::vector<double> &pixdim)
    {
        Affine result = identity();
        for (std::size_t i=0; i<3 && i<pixdim.size(); i++)
            result(static_cast<int>(i), static_cast<int>(i)) = pixdim[i];
        return result;
    }

    double & operator() (const int i, const int j) { return values[i*4 + j]; }
    double operator() (const int i, const int j) const { return values[i*4 + j]; }

    const std::array<double,16> & data () const { return values; }

    // Apply to a position, which is treated as having an implicit fourth
    // component of one so that the translation is included
    Point multiply (const Point &p) const
    {
        Point result;
        for (int i=0; i<3; i++)
            result[i] = (*this)(i,0)*p[0] + (*this)(i,1)*p[1] + (*this)(i,2)*p[2] + (*this)(i,3);
        return result;
    }

    bool isAffine () const
    {
        const double tolerance = 1e-10;
        return std::fabs((*this)(3,0)) < tolerance && std::fabs((*this)(3,1)) < tolerance
            && std::fabs((*this)(3,2)) < tolerance && std::fabs((*this)(3,3) - 1.0) < tolerance;
    }

    double determinant3 () const;

    // Inverse of [R t; 0 1] is [R^-1, -R^-1 t; 0 1], so only the 3x3 block has
    // to be inverted. Throws if that block is singular
    Affine inverse () const;
};

// The geometry of the space an image is embedded within. Deliberately free of
// any dependency on a file format: it holds only what the mapping needs, and
// NIfTI or other interop is layered on top
class ImageSpace
{
public:
    int spatial;
    std::vector<double> pixdim;
    Affine xform;
    std::string spaceUnit, timeUnit;

    ImageSpace ()
        : spatial(0), xform(Affine::identity()), spaceUnit("unknown"), timeUnit("unknown") {}

    ImageSpace (const int spatial, const std::vector<double> &pixdim)
        : spatial(spatial), pixdim(pixdim), xform(Affine::scaling(pixdim)),
          spaceUnit("unknown"), timeUnit("unknown") {}

    ImageSpace (const int spatial, const std::vector<double> &pixdim, const Affine &xform)
        : spatial(spatial), pixdim(pixdim), xform(xform),
          spaceUnit("unknown"), timeUnit("unknown") {}

    // Convert a point of the given type to voxel coordinates
    Point toVoxel (const Point &p, const PointType type) const;

    // The reverse: voxel coordinates to a point of the given type
    Point fromVoxel (const Point &p, const PointType type) const;
};

// Rounding is kept separate from coordinate conversion. Bounds are optional
// and only consulted by the probabilistic strategy, to avoid selecting a
// location off the end of the image
//
// The generator is passed in rather than being global, so this is safe to call
// from a worker thread provided each thread owns its generator. It is only
// consulted by the probabilistic strategy, and may be null otherwise
Point roundLocation (const Point &p, const RoundingType round,
                     const std::vector<std::size_t> *bounds = nullptr,
                     RandomGenerator *generator = nullptr);

// Definitions are inline and live here rather than in a source file, so
// that a package linking to imply needs only the headers
inline double Affine::determinant3 () const
{
    const Affine &m = *this;
    return m(0,0) * (m(1,1)*m(2,2) - m(1,2)*m(2,1))
         - m(0,1) * (m(1,0)*m(2,2) - m(1,2)*m(2,0))
         + m(0,2) * (m(1,0)*m(2,1) - m(1,1)*m(2,0));
}

inline Affine Affine::inverse () const
{
    if (!isAffine())
        Rcpp::stop("Only affine transforms can be inverted this way");

    const Affine &m = *this;
    const double det = determinant3();
    if (std::fabs(det) < 1e-12)
        Rcpp::stop("Transform matrix is singular and cannot be inverted");

    Affine result;

    // Inverse of the 3x3 block, by the adjugate
    result(0,0) = (m(1,1)*m(2,2) - m(1,2)*m(2,1)) / det;
    result(0,1) = (m(0,2)*m(2,1) - m(0,1)*m(2,2)) / det;
    result(0,2) = (m(0,1)*m(1,2) - m(0,2)*m(1,1)) / det;
    result(1,0) = (m(1,2)*m(2,0) - m(1,0)*m(2,2)) / det;
    result(1,1) = (m(0,0)*m(2,2) - m(0,2)*m(2,0)) / det;
    result(1,2) = (m(0,2)*m(1,0) - m(0,0)*m(1,2)) / det;
    result(2,0) = (m(1,0)*m(2,1) - m(1,1)*m(2,0)) / det;
    result(2,1) = (m(0,1)*m(2,0) - m(0,0)*m(2,1)) / det;
    result(2,2) = (m(0,0)*m(1,1) - m(0,1)*m(1,0)) / det;

    // Translation becomes -R^-1 t
    for (int i=0; i<3; i++)
        result(i,3) = -(result(i,0)*m(0,3) + result(i,1)*m(1,3) + result(i,2)*m(2,3));

    result(3,3) = 1.0;
    return result;
}

inline Point ImageSpace::toVoxel (const Point &p, const PointType type) const
{
    Point result = p;

    switch (type)
    {
        case PointType::voxel:
        break;

        case PointType::scaled:
        for (int i=0; i<3; i++)
        {
            const double scale = (i < static_cast<int>(pixdim.size()) ? std::fabs(pixdim[i]) : 1.0);
            result[i] = (scale > 0.0 ? p[i] / scale : p[i]);
        }
        break;

        // The stored transform maps voxel coordinates to world coordinates, so
        // going the other way needs its inverse. The version this was ported
        // from applied the forward transform here, which was a bug
        case PointType::world:
        result = xform.inverse().multiply(p);
        break;
    }

    return result;
}

inline Point ImageSpace::fromVoxel (const Point &p, const PointType type) const
{
    Point result = p;

    switch (type)
    {
        case PointType::voxel:
        break;

        case PointType::scaled:
        for (int i=0; i<3; i++)
        {
            const double scale = (i < static_cast<int>(pixdim.size()) ? std::fabs(pixdim[i]) : 1.0);
            result[i] = p[i] * scale;
        }
        break;

        case PointType::world:
        result = xform.multiply(p);
        break;
    }

    return result;
}

inline Point roundLocation (const Point &p, const RoundingType round, const std::vector<std::size_t> *bounds,
                     RandomGenerator *generator)
{
    Point result = p;

    switch (round)
    {
        case RoundingType::none:
        break;

        // nearbyint rather than round, so that a coordinate falling exactly
        // halfway breaks to even, matching R's round()
        case RoundingType::conventional:
        for (int i=0; i<3; i++)
            result[i] = std::nearbyint(p[i]);
        break;

        case RoundingType::probabilistic:
        if (generator == nullptr)
            Rcpp::stop("Probabilistic rounding requires a random number generator");

        for (int i=0; i<3; i++)
        {
            const double ceiling = std::ceil(p[i]);
            const double floor = std::floor(p[i]);
            const double distance = p[i] - floor;

            // Sample in proportion to proximity, unless that would step off the
            // end of the image
            const double sample = generator->uniform();
            const bool beyondEnd = (bounds != nullptr && i < static_cast<int>(bounds->size())
                                    && ceiling >= static_cast<double>((*bounds)[i]));
            const bool chooseFloor = (sample > distance && floor >= 0.0) || beyondEnd;
            result[i] = chooseFloor ? floor : ceiling;
        }
        break;
    }

    return result;
}

} // namespace imply

#endif
