#ifndef SWE_TRUE_SOLUTION_FUNCTIONS_HPP
#define SWE_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> true_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];
    Utilities::ignore(t);
  
    constexpr double alpha = 1.e-7;
    constexpr double X = 1;
    constexpr double Y = -0.41884;
    const double w = sqrt(8*Global::g*alpha);
    const double tau = 2*PI/w;

    double r2 = x*x + y*y;
    double bathymetry = alpha * r2;
    
    double water_column_height = 1/(X+Y*cos(w*t)) + alpha*(Y*Y-X*X)*r2/((X+Y*cos(w*t))*(X+Y*cos(w*t)));
    double true_ze = water_column_height - bathymetry;
    double true_qx = -Y*w*sin(w*t)*x/(2*(X+Y*cos(w*t))) * water_column_height;
    double true_qy = -Y*w*sin(w*t)*y/(2*(X+Y*cos(w*t))) * water_column_height;
  
    StatVector<double, SWE::n_variables> true_q{true_ze, true_qx, true_qy};

    return true_q;
}
}

#endif
