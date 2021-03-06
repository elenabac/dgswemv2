#ifndef EHDG_SWE_IS_INTERNAL_HPP
#define EHDG_SWE_IS_INTERNAL_HPP

namespace SWE {
namespace EHDG {
namespace ISP {
class Internal {
  public:
    template <typename InterfaceType>
    void Initialize(InterfaceType& intface) {} /*nothing to initialize*/

    template <typename EdgeInterfaceType>
    void ComputeGlobalKernels(EdgeInterfaceType& edge_int);

    template <typename EdgeInterfaceType>
    void ComputeNumericalFlux(EdgeInterfaceType& edge_int);
};

template <typename EdgeInterfaceType>
void Internal::ComputeGlobalKernels(EdgeInterfaceType& edge_int) {
    auto& edge_internal = edge_int.edge_data.edge_internal;

    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    uint gp_ex;
    for (uint gp = 0; gp < edge_int.edge_data.get_ngp(); ++gp) {
        gp_ex = edge_int.edge_data.get_ngp() - gp - 1;

        column(edge_internal.rhs_global_kernel_at_gp, gp) = column(boundary_in.Fn_at_gp, gp);
        column(edge_internal.rhs_global_kernel_at_gp, gp) += column(boundary_ex.Fn_at_gp, gp_ex);
    }

    // Add tau terms
    add_kernel_tau_terms_intface_LF(edge_int);
}

template <typename EdgeInterfaceType>
void Internal::ComputeNumericalFlux(EdgeInterfaceType& edge_int) {
    auto& boundary_in = edge_int.interface.data_in.boundary[edge_int.interface.bound_id_in];
    auto& boundary_ex = edge_int.interface.data_ex.boundary[edge_int.interface.bound_id_ex];

    boundary_in.F_hat_at_gp = boundary_in.Fn_at_gp;
    boundary_ex.F_hat_at_gp = boundary_ex.Fn_at_gp;

    // Add tau terms
    add_F_hat_tau_terms_intface_LF(edge_int);
}
}
}
}

#endif