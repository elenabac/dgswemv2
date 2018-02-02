#ifndef SHAPES_2D_HPP
#define SHAPES_2D_HPP

#include "general_definitions.hpp"

namespace Shape {
class StraightTriangle : public Shape<2> {
  public:
    StraightTriangle(const std::vector<Point<2>>& nodal_coordinates);

    std::vector<uint> GetBoundaryNodeID(const uint bound_id, const std::vector<uint> node_ID);

    Point<2> GetBarycentricCoordinates() const;
    std::vector<Point<2>> GetMidpointCoordinates() const;

    std::vector<double> GetJdet(const std::vector<Point<2>>& points) const;
    Array3D<double> GetJinv(const std::vector<Point<2>>& points) const;
    std::vector<double> GetSurfaceJ(const uint bound_id, const std::vector<Point<2>>& points) const;
    Array2D<double> GetSurfaceNormal(const uint bound_id, const std::vector<Point<2>>& points) const;

    Array2D<double> GetPsi(const std::vector<Point<2>>& points);
    Array3D<double> GetDPsi(const std::vector<Point<2>>& points);

    Array2D<double> GetBoundaryPsi(const uint bound_id, const std::vector<Point<1>>& points);

    std::vector<Point<2>> LocalToGlobalCoordinates(const std::vector<Point<2>>& points);

    void GetVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) const;

  private:
    std::vector<Point<2>> nodal_coordinates;

#ifdef HAS_HPX
  public:
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

#ifdef HAS_HPX
template <typename Archive>
void StraightTriangle::serialize(Archive& ar, unsigned) {
    ar & nodal_coordinates;
}
#endif
}

#endif
