#ifndef _IMPLY_SPACE_H_
#define _IMPLY_SPACE_H_

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
class randomGenerator
{
protected:
    std::mt19937_64 engine;
    std::uniform_real_distribution<double> distribution;

public:
    explicit randomGenerator (const std::uint64_t seed)
        : engine(seed), distribution(0.0, 1.0) {}

    // Uniform on [0,1), matching the convention of R's unif_rand()
    double uniform () { return distribution(engine); }

    void reseed (const std::uint64_t seed) { engine.seed(seed); }
};

// Location conventions: voxel-indexed, scaled for voxel dimensions only (as
// with a diagonal xform), or world coordinates fully respecting the xform
enum class pointType { voxel, scaled, world };

// Rounding strategies: none, standard for nearest-neighbour, or probabilistic
// for stochastic nearest neighbour (probabilities proportional to distance)
enum class roundingType { none, conventional, probabilistic };

typedef std::array<double,3> point;

// A 4x4 affine transform, stored row-major. Only the affine case is supported,
// meaning the final row is implicitly (0,0,0,1), which is what makes the
// inverse cheap and exact
class affine
{
protected:
    std::array<double,16> values;

public:
    affine () { values.fill(0.0); }

    static affine identity ()
    {
        affine result;
        for (int i=0; i<4; i++)
            result(i,i) = 1.0;
        return result;
    }

    // A diagonal transform built from voxel dimensions, used when an image
    // carries no explicit transform of its own
    static affine scaling (const std::vector<double> &pixdim)
    {
        affine result = identity();
        for (std::size_t i=0; i<3 && i<pixdim.size(); i++)
            result(static_cast<int>(i), static_cast<int>(i)) = pixdim[i];
        return result;
    }

    double & operator() (const int i, const int j) { return values[i*4 + j]; }
    double operator() (const int i, const int j) const { return values[i*4 + j]; }

    const std::array<double,16> & data () const { return values; }

    // Apply to a position, which is treated as having an implicit fourth
    // component of one so that the translation is included
    point multiply (const point &p) const
    {
        point result;
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
    affine inverse () const;
};

// The geometry of the space an image is embedded within. Deliberately free of
// any dependency on a file format: it holds only what the mapping needs, and
// NIfTI or other interop is layered on top
class imageSpace
{
public:
    int spatial;
    std::vector<double> pixdim;
    affine xform;
    std::string spaceUnit, timeUnit;

    imageSpace ()
        : spatial(0), xform(affine::identity()), spaceUnit("unknown"), timeUnit("unknown") {}

    imageSpace (const int spatial, const std::vector<double> &pixdim)
        : spatial(spatial), pixdim(pixdim), xform(affine::scaling(pixdim)),
          spaceUnit("unknown"), timeUnit("unknown") {}

    imageSpace (const int spatial, const std::vector<double> &pixdim, const affine &xform)
        : spatial(spatial), pixdim(pixdim), xform(xform),
          spaceUnit("unknown"), timeUnit("unknown") {}

    // Three-letter code naming, for each voxel axis, the anatomical direction
    // in which its index increases. Derived from the transform matrix directly,
    // so no quaternion representation is needed
    std::string orientation () const;

    // Convert a point of the given type to voxel coordinates
    point toVoxel (const point &p, const pointType type) const;

    // The reverse: voxel coordinates to a point of the given type
    point fromVoxel (const point &p, const pointType type) const;
};

// Rounding is kept separate from coordinate conversion, which the original
// conflated. Bounds are optional and only consulted by the probabilistic
// strategy, to avoid selecting a location off the end of the image.
//
// The generator is passed in rather than being global, so this is safe to call
// from a worker thread provided each thread owns its generator. It is only
// consulted by the probabilistic strategy, and may be null otherwise
point roundLocation (const point &p, const roundingType round,
                     const std::vector<std::size_t> *bounds = nullptr,
                     randomGenerator *generator = nullptr);

} // namespace imply

#endif
