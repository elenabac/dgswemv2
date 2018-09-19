#ifndef IHDG_SWE_PROC_OMPI_STAGE_HPP
#define IHDG_SWE_PROC_OMPI_STAGE_HPP

#include "general_definitions.hpp"

#include "ihdg_swe_kernels_processor.hpp"
#include "ihdg_swe_proc_ompi_sol_glob_prob.hpp"

namespace SWE {
namespace IHDG {
template <typename OMPISimUnitType>
void Problem::stage_ompi(std::vector<std::unique_ptr<OMPISimUnitType>>& sim_units) {
#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            const uint stage = sim_units[su_id]->stepper.GetStage();

            auto& state      = elt.data.state[stage + 1];
            auto& state_prev = elt.data.state[stage];
            auto& internal   = elt.data.internal;

            state.q = state_prev.q;

            internal.q_prev_at_gp = elt.ComputeUgp(state_prev.q);
        });
    }

    uint iter = 0;
    while (true) {
        iter++;

#pragma omp parallel for
        for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
            /* Local Step */
            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::local_volume_kernel(sim_units[su_id]->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachElement(
                [&sim_units, su_id](auto& elt) { Problem::local_source_kernel(sim_units[su_id]->stepper, elt); });

            sim_units[su_id]->discretization.mesh.CallForEachInterface([&sim_units, su_id](auto& intface) {
                Problem::local_interface_kernel(sim_units[su_id]->stepper, intface);
            });

            sim_units[su_id]->discretization.mesh.CallForEachBoundary(
                [&sim_units, su_id](auto& bound) { Problem::local_boundary_kernel(sim_units[su_id]->stepper, bound); });

            sim_units[su_id]->discretization.mesh.CallForEachDistributedBoundary([&sim_units, su_id](auto& dbound) {
                Problem::local_distributed_boundary_kernel(sim_units[su_id]->stepper, dbound);
            });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&sim_units, su_id](auto& edge_int) {
                    Problem::local_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&sim_units, su_id](auto& edge_bound) {
                    Problem::local_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&sim_units, su_id](auto& edge_dbound) {
                    Problem::local_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
                });
            /* Local Step */

            /* Global Step */
            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeInterface(
                [&sim_units, su_id](auto& edge_int) {
                    Problem::global_edge_interface_kernel(sim_units[su_id]->stepper, edge_int);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeBoundary(
                [&sim_units, su_id](auto& edge_bound) {
                    Problem::global_edge_boundary_kernel(sim_units[su_id]->stepper, edge_bound);
                });

            sim_units[su_id]->discretization.mesh_skeleton.CallForEachEdgeDistributed(
                [&sim_units, su_id](auto& edge_dbound) {
                    Problem::global_edge_distributed_kernel(sim_units[su_id]->stepper, edge_dbound);
                });
            /* Global Step */
        }

        bool converged = ompi_solve_global_problem(sim_units);

        if (converged) {
            break;
        }

        if (iter == 100) {
            break;
        }
    }

#pragma omp parallel for
    for (uint su_id = 0; su_id < sim_units.size(); ++su_id) {
        sim_units[su_id]->discretization.mesh.CallForEachElement([&sim_units, su_id](auto& elt) {
            bool nan_found = Problem::scrutinize_solution_kernel(sim_units[su_id]->stepper, elt);

            if (nan_found)
                MPI_Abort(MPI_COMM_WORLD, 0);
        });

        ++(sim_units[su_id]->stepper);
    }
}
}
}

#endif