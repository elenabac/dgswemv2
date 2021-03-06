#ifndef RKDG_SWE_DATA_BOUNDARY_HPP
#define RKDG_SWE_DATA_BOUNDARY_HPP

namespace SWE {
namespace RKDG {
struct Boundary {
    Boundary() = default;
    Boundary(const uint ngp)
        : q_at_gp(SWE::n_variables, ngp), aux_at_gp(SWE::n_auxiliaries, ngp), F_hat_at_gp(SWE::n_variables, ngp) {}

    HybMatrix<double, SWE::n_variables> q_at_gp;
    HybMatrix<double, SWE::n_auxiliaries> aux_at_gp;

    HybMatrix<double, SWE::n_variables> F_hat_at_gp;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & q_at_gp
            & aux_at_gp
            & F_hat_at_gp;
        // clang-format on
    }
#endif
};
}
}

#endif
