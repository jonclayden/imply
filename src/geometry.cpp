#include <Rcpp.h>

#include <memory>

#include "imply/Space.h"

using namespace imply;

namespace {

Affine affineFrom (const Rcpp::NumericMatrix &m)
{
    if (m.nrow() != 4 || m.ncol() != 4)
        Rcpp::stop("Transform matrix must be 4x4, not %dx%d", m.nrow(), m.ncol());

    Affine result;
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            result(i,j) = m(i,j);
    }
    return result;
}

Rcpp::NumericMatrix affineTo (const Affine &a)
{
    Rcpp::NumericMatrix result(4, 4);
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            result(i,j) = a(i,j);
    }
    return result;
}

PointType parsePointType (const std::string &name)
{
    if (name == "voxel") return PointType::voxel;
    if (name == "scaled") return PointType::scaled;
    if (name == "world") return PointType::world;
    Rcpp::stop("Point type should be \"voxel\", \"scaled\" or \"world\", not \"%s\"", name);
}

RoundingType parseRoundingType (const std::string &name)
{
    if (name == "none") return RoundingType::none;
    if (name == "conventional") return RoundingType::conventional;
    if (name == "probabilistic") return RoundingType::probabilistic;
    Rcpp::stop("Rounding type should be \"none\", \"conventional\" or \"probabilistic\", not \"%s\"", name);
}

ImageSpace spaceFrom (const Rcpp::NumericMatrix &xform, const Rcpp::NumericVector &pixdim)
{
    const std::vector<double> dims(pixdim.begin(), pixdim.end());
    return ImageSpace(static_cast<int>(dims.size()), dims, affineFrom(xform));
}

// Points arrive as a matrix with one row per point and three columns
Rcpp::NumericMatrix convertPoints (const Rcpp::NumericMatrix &locs, const ImageSpace &space,
                                   const bool toVoxel, const PointType type)
{
    if (locs.ncol() != 3)
        Rcpp::stop("Point matrix must have three columns, not %d", locs.ncol());

    Rcpp::NumericMatrix result(locs.nrow(), 3);
    for (R_xlen_t i=0; i<locs.nrow(); i++)
    {
        const Point source = { locs(i,0), locs(i,1), locs(i,2) };
        const Point converted = (toVoxel ? space.toVoxel(source, type) : space.fromVoxel(source, type));
        for (int j=0; j<3; j++)
            result(i,j) = converted[j];
    }

    return result;
}

} // anonymous namespace

// [[Rcpp::export]]
std::string orientationFromXform (Rcpp::NumericMatrix xform)
{
    ImageSpace space;
    space.xform = affineFrom(xform);
    return space.orientation();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix invertXform (Rcpp::NumericMatrix xform)
{
    return affineTo(affineFrom(xform).inverse());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pointsToVoxel (Rcpp::NumericMatrix locs, Rcpp::NumericMatrix xform,
                                   Rcpp::NumericVector pixdim, std::string type = "world")
{
    return convertPoints(locs, spaceFrom(xform, pixdim), true, parsePointType(type));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pointsFromVoxel (Rcpp::NumericMatrix locs, Rcpp::NumericMatrix xform,
                                     Rcpp::NumericVector pixdim, std::string type = "world")
{
    return convertPoints(locs, spaceFrom(xform, pixdim), false, parsePointType(type));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix roundPoints (Rcpp::NumericMatrix locs, std::string round = "conventional",
                                 Rcpp::Nullable<Rcpp::NumericVector> bounds = R_NilValue)
{
    if (locs.ncol() != 3)
        Rcpp::stop("Point matrix must have three columns, not %d", locs.ncol());

    const RoundingType strategy = parseRoundingType(round);

    std::vector<std::size_t> extents;
    if (bounds.isNotNull())
    {
        const Rcpp::NumericVector values(bounds.get());
        extents.assign(values.begin(), values.end());
    }
    const std::vector<std::size_t> *extentsPtr = (extents.empty() ? nullptr : &extents);

    // The generator is our own, so nothing here touches R's global RNG beyond
    // drawing a single seed. That keeps set.seed() in charge of reproducibility
    // while leaving the rounding itself safe to call from a worker thread
    std::unique_ptr<RandomGenerator> generator;
    if (strategy == RoundingType::probabilistic)
    {
        GetRNGstate();
        const std::uint64_t seed = static_cast<std::uint64_t>(unif_rand() * 9007199254740992.0);
        PutRNGstate();
        generator.reset(new RandomGenerator(seed));
    }

    Rcpp::NumericMatrix result(locs.nrow(), 3);
    for (R_xlen_t i=0; i<locs.nrow(); i++)
    {
        const Point source = { locs(i,0), locs(i,1), locs(i,2) };
        const Point rounded = roundLocation(source, strategy, extentsPtr, generator.get());
        for (int j=0; j<3; j++)
            result(i,j) = rounded[j];
    }

    return result;
}
