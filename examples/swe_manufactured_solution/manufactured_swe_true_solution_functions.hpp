#ifndef SWE_TRUE_SOLUTION_FUNCTIONS_HPP
#define SWE_TRUE_SOLUTION_FUNCTIONS_HPP

#include "utilities/ignore.hpp"

namespace SWE {
inline StatVector<double, SWE::n_variables> true_q(const double t, const Point<2>& pt) {
    const double x = pt[GlobalCoord::x];
    const double y = pt[GlobalCoord::y];

    /*
    constexpr double x1 = -PI;
    constexpr double x2 = PI;
    constexpr double y1 = -PI;
    constexpr double y2 = PI;
    */
    constexpr double x1 = 40000.;
    constexpr double x2 = 83200.;
    constexpr double y1 = 10000.;
    constexpr double y2 = 53200.;
    
    constexpr double Ho = 2.;
    constexpr double zo = 0.25;
    Utilities::ignore(Ho);

    constexpr double tau = 3456.;
    constexpr double w   = 2*PI/43200.;

    double true_ze = 2 * zo * cos(w * (x - x1)) * cos(w * (y - y1)) *
                     cos(w * (t + tau)) / (cos(w * (x2 - x1)) * cos(w * (y2 - y1))) + Ho;

    double true_qx = zo * sin(w * (x - x1)) * cos(w * (y - y1)) * sin(w * (t + tau)) /
                     (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    double true_qy = zo * cos(w * (x - x1)) * sin(w * (y - y1)) * sin(w * (t + tau)) /
                     (cos(w * (x2 - x1)) * cos(w * (y2 - y1)));

    /*
    double true_ze = exp(sin(3 * x) * sin(3 * y) - sin(3 * t));
    double true_qx = cos(x - 4 * t);
    double true_qy = sin(y + 4 * t);
    */
    
    StatVector<double, SWE::n_variables> true_q{true_ze, true_qx, true_qy};

    return true_q;
}
}

#endif
