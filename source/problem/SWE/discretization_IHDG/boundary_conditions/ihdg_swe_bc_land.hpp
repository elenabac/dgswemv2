#ifndef IHDG_SWE_BC_LAND_HPP
#define IHDG_SWE_BC_LAND_HPP

namespace SWE {
namespace IHDG {
namespace BC {
class Land {
  public:
    template <typename BoundaryType>
    void Initialize(BoundaryType& bound) {} /*nothing to initialize*/

    template <typename StepperType, typename EdgeBoundaryType>
    void ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound);
};

template <typename StepperType, typename EdgeBoundaryType>
void Land::ComputeGlobalKernels(const StepperType& stepper, EdgeBoundaryType& edge_bound) {
    auto& edge_internal = edge_bound.edge_data.edge_internal;

    auto& boundary = edge_bound.boundary.data.boundary[edge_bound.boundary.bound_id];

    double qn;
    double nx, ny;

    StatVector<double, SWE::n_variables* SWE::n_variables> I_vector = IdentityVector<double>(SWE::n_variables);

    for (uint gp = 0; gp < edge_bound.edge_data.get_ngp(); ++gp) {
        nx = edge_bound.boundary.surface_normal(GlobalCoord::x, gp);
        ny = edge_bound.boundary.surface_normal(GlobalCoord::y, gp);

        qn = boundary.q_at_gp(SWE::Variables::qx, gp) * nx + boundary.q_at_gp(SWE::Variables::qy, gp) * ny;

        column(boundary.delta_global_kernel_at_gp, gp) = -I_vector;

        boundary.delta_global_kernel_at_gp(JacobianVariables::qx_qx, gp) += nx * nx;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qx_qy, gp) += nx * ny;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qy_qx, gp) += nx * ny;
        boundary.delta_global_kernel_at_gp(JacobianVariables::qy_qy, gp) += ny * ny;

        column(edge_internal.delta_hat_global_kernel_at_gp, gp) = I_vector;

        column(edge_internal.rhs_global_kernel_at_gp, gp) =
            column(edge_internal.q_hat_at_gp, gp) - column(boundary.q_at_gp, gp);
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qx, gp) += qn * nx;
        edge_internal.rhs_global_kernel_at_gp(SWE::Variables::qy, gp) += qn * ny;
    }
}
}
}
}

#endif