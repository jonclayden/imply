#include <Rcpp.h>

#include "Space.h"

namespace imply {

double affine::determinant3 () const
{
    const affine &m = *this;
    return m(0,0) * (m(1,1)*m(2,2) - m(1,2)*m(2,1))
         - m(0,1) * (m(1,0)*m(2,2) - m(1,2)*m(2,0))
         + m(0,2) * (m(1,0)*m(2,1) - m(1,1)*m(2,0));
}

affine affine::inverse () const
{
    if (!isAffine())
        Rcpp::stop("Only affine transforms can be inverted this way");

    const affine &m = *this;
    const double det = determinant3();
    if (std::fabs(det) < 1e-12)
        Rcpp::stop("Transform matrix is singular and cannot be inverted");

    affine result;

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

std::string imageSpace::orientation () const
{
    // Column j of the 3x3 block is the world-space direction along which voxel
    // axis j increases, so each voxel axis has to be matched to the anatomical
    // axis it aligns with most strongly.
    //
    // Columns are first normalised to unit length, so that anisotropic voxel
    // dimensions cannot outweigh direction, and the assignment is then chosen
    // by exhaustive search over all six permutations. A greedy nearest-axis
    // assignment is not equivalent: it disagrees with the NIfTI reference
    // implementation on roughly one oblique transform in fifteen
    static const char codes[3][2] = { {'L','R'}, {'P','A'}, {'I','S'} };
    static const int permutations[6][3] = { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} };

    double q[3][3];
    for (int j=0; j<3; j++)
    {
        double norm = 0.0;
        for (int i=0; i<3; i++)
            norm += xform(i,j) * xform(i,j);
        norm = std::sqrt(norm);

        for (int i=0; i<3; i++)
            q[i][j] = (norm > 0.0 ? xform(i,j) / norm : xform(i,j));
    }

    int best = -1;
    double bestScore = -1.0;
    for (int p=0; p<6; p++)
    {
        double score = 0.0;
        for (int j=0; j<3; j++)
            score += std::fabs(q[permutations[p][j]][j]);

        if (score > bestScore)
        {
            bestScore = score;
            best = p;
        }
    }

    std::string result(3, '?');
    if (best < 0)
        return result;

    for (int j=0; j<3; j++)
    {
        const int i = permutations[best][j];
        result[j] = codes[i][q[i][j] > 0.0 ? 1 : 0];
    }

    return result;
}

point imageSpace::toVoxel (const point &p, const pointType type) const
{
    point result = p;

    switch (type)
    {
        case pointType::voxel:
        break;

        case pointType::scaled:
        for (int i=0; i<3; i++)
        {
            const double scale = (i < static_cast<int>(pixdim.size()) ? std::fabs(pixdim[i]) : 1.0);
            result[i] = (scale > 0.0 ? p[i] / scale : p[i]);
        }
        break;

        // The stored transform maps voxel coordinates to world coordinates, so
        // going the other way needs its inverse. The version this was ported
        // from applied the forward transform here, which was a bug
        case pointType::world:
        result = xform.inverse().multiply(p);
        break;
    }

    return result;
}

point imageSpace::fromVoxel (const point &p, const pointType type) const
{
    point result = p;

    switch (type)
    {
        case pointType::voxel:
        break;

        case pointType::scaled:
        for (int i=0; i<3; i++)
        {
            const double scale = (i < static_cast<int>(pixdim.size()) ? std::fabs(pixdim[i]) : 1.0);
            result[i] = p[i] * scale;
        }
        break;

        case pointType::world:
        result = xform.multiply(p);
        break;
    }

    return result;
}

point roundLocation (const point &p, const roundingType round, const std::vector<std::size_t> *bounds,
                     randomGenerator *generator)
{
    point result = p;

    switch (round)
    {
        case roundingType::none:
        break;

        // nearbyint rather than round, so that a coordinate falling exactly
        // halfway breaks to even, matching R's round(). The code this was
        // ported from used std::round, which breaks away from zero and so
        // disagrees with R on values like 4.5
        case roundingType::conventional:
        for (int i=0; i<3; i++)
            result[i] = std::nearbyint(p[i]);
        break;

        case roundingType::probabilistic:
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
