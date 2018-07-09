#ifndef RKDG_SWE_DATA_BOUNDARY_HPP
#define RKDG_SWE_DATA_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace SWE {
namespace RKDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp) : q_at_gp(ngp), bath_at_gp(ngp), h_at_gp(ngp), F_hat_at_gp(ngp) {}

    std::vector<Vector<double, SWE::n_variables>> q_at_gp;
    std::vector<double> bath_at_gp;
    std::vector<double> h_at_gp;

    std::vector<Vector<double, SWE::n_variables>> F_hat_at_gp;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        /*ar  & q_at_gp
            & bath_at_gp
            & h_at_gp
            & F_hat_at_gp;*/
        // clang-format on
    }
#endif
};
}
}

#endif
