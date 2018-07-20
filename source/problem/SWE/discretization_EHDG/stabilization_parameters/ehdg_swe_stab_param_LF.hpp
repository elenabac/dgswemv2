#ifndef EHDG_SWE_STAB_PARAM_LF_HPP
#define EHDG_SWE_STAB_PARAM_LF_HPP

#include "problem/SWE/swe_definitions.hpp"

namespace SWE {
namespace EHDG {
template <typename EdgeInterfaceType>
inline void add_kernel_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    StatVector<double, SWE::n_variables> del_q;
    StatVector<double, SWE::n_variables> dtau_dq_hat;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dq_hat[SWE::Variables::ze] =
            std::sqrt(Global::g / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]) / 2.0 -
            sgn * un_hat / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dq_hat[SWE::Variables::qx] = sgn * nx / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dq_hat[SWE::Variables::qy] = sgn * ny / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        del_q = boundary_in.q_at_gp[gp] + boundary_ex.q_at_gp[gp_ex] - 2.0 * edge_internal.q_hat_at_gp[gp];

        /* delta kernels */

        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::ze_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::ze];
        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::qx_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::qx];
        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::qy_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::qy];

        edge_internal.delta_hat_global_kernel_at_gp[gp] += -2 * tau * I_vector;

        /* RHS kernels */

        edge_internal.rhs_global_kernel_at_gp[gp] += tau * del_q;
    }
}

template <typename EdgeDistributedType>
inline void add_kernel_tau_terms_dbound_LF(EdgeDistributedType& edge_dbound) {
    auto& edge_internal = edge_dbound.edge_data.edge_internal;

    auto& boundary = edge_dbound.boundary.data.boundary[edge_dbound.boundary.bound_id];

    double u_hat, v_hat, un_hat;
    double nx, ny;
    double tau, dtau_dze_hat, dtau_dqx_hat, dtau_dqy_hat;
    double sgn;

    StatVector<double, SWE::n_variables> del_q;
    StatVector<double, SWE::n_variables> dtau_dq_hat;

    StatVector<double, SWE::n_variables> q_ex;
    StatVector<double, SWE::n_variables> Fn_ex;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < edge_dbound.edge_data.get_ngp(); ++gp) {
        edge_dbound.boundary.boundary_condition.exchanger.GetEX(gp, q_ex, Fn_ex);

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_dbound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_dbound.boundary.surface_normal[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        if (0.0 < un_hat) {
            sgn = 1.0;
        } else if (un_hat < 0.0) {
            sgn = -1.0;
        } else {
            sgn = 0.0;
        }

        dtau_dq_hat[SWE::Variables::ze] =
            std::sqrt(Global::g / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]) / 2.0 -
            sgn * un_hat / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dq_hat[SWE::Variables::qx] = sgn * nx / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        dtau_dq_hat[SWE::Variables::qy] = sgn * ny / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        del_q = boundary.q_at_gp[gp] + q_ex - 2.0 * edge_internal.q_hat_at_gp[gp];

        /* delta kernels */

        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::ze_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::ze];
        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::qx_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::qx];
        subvector(edge_internal.delta_hat_global_kernel_at_gp[gp], SWE::JacobianVariables::qy_ze, SWE::n_variables) =
            dtau_dq_hat * del_q[SWE::Variables::qy];

        edge_internal.delta_hat_global_kernel_at_gp[gp] += -2 * tau * I_vector;

        /* RHS kernels */

        edge_internal.rhs_global_kernel_at_gp[gp] += tau * del_q;
    }
}

template <typename EdgeInterfaceType>
inline void add_F_hat_tau_terms_intface_LF(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_int.interface.surface_normal_in[gp][GlobalCoord::x];
        ny = edge_int.interface.surface_normal_in[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        boundary_in.F_hat_at_gp[gp] += tau * (boundary_in.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp]);
        boundary_ex.F_hat_at_gp[gp_ex] += tau * (boundary_ex.q_at_gp[gp_ex] - edge_internal.q_hat_at_gp[gp]);
    }
}

template <typename EdgeBoundaryType>
inline void add_F_hat_tau_terms_bound_LF(EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double nx, ny;
    double u_hat, v_hat, un_hat;
    double tau;

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        u_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qx] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];
        v_hat =
            edge_internal.q_hat_at_gp[gp][SWE::Variables::qy] / edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h];

        nx = edge_bound.boundary.surface_normal[gp][GlobalCoord::x];
        ny = edge_bound.boundary.surface_normal[gp][GlobalCoord::y];

        un_hat = u_hat * nx + v_hat * ny;

        tau = std::abs(un_hat) + std::sqrt(Global::g * edge_internal.aux_hat_at_gp[gp][SWE::Auxiliaries::h]);

        boundary.F_hat_at_gp[gp] += tau * (boundary.q_at_gp[gp] - edge_internal.q_hat_at_gp[gp]);
    }
}
}
}

#endif