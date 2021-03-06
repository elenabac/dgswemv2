#ifndef RKDG_SWE_PROC_HPX_STAGE_HPP
#define RKDG_SWE_PROC_HPX_STAGE_HPP

#include "rkdg_swe_kernels_processor.hpp"
#include "problem/SWE/problem_slope_limiter/swe_CS_sl_hpx.hpp"

namespace SWE {
namespace RKDG {
template <typename HPXSimUnitType>
auto Problem::stage_hpx(HPXSimUnitType* sim_unit) {
    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Current (time, stage): (" << sim_unit->stepper.GetTimeAtCurrentStage() << ','
                                      << sim_unit->stepper.GetStage() << ')' << std::endl;

        sim_unit->writer.GetLogFile() << "Exchanging data" << std::endl;
    }

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(CommTypes::bound_state, sim_unit->stepper.GetTimestamp());

    sim_unit->discretization.mesh.CallForEachDistributedBoundary(
        [sim_unit](auto& dbound) { Problem::distributed_boundary_send_kernel(sim_unit->stepper, dbound); });

    sim_unit->communicator.SendAll(CommTypes::bound_state, sim_unit->stepper.GetTimestamp());

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    sim_unit->discretization.mesh.CallForEachElement(
        [sim_unit](auto& elt) { Problem::volume_kernel(sim_unit->stepper, elt); });

    sim_unit->discretization.mesh.CallForEachElement(
        [sim_unit](auto& elt) { Problem::source_kernel(sim_unit->stepper, elt); });

    sim_unit->discretization.mesh.CallForEachInterface(
        [sim_unit](auto& intface) { Problem::interface_kernel(sim_unit->stepper, intface); });

    sim_unit->discretization.mesh.CallForEachBoundary(
        [sim_unit](auto& bound) { Problem::boundary_kernel(sim_unit->stepper, bound); });

    if (sim_unit->writer.WritingVerboseLog()) {
        sim_unit->writer.GetLogFile() << "Finished work before receive" << std::endl
                                      << "Starting to wait on receive with timestamp: "
                                      << sim_unit->stepper.GetTimestamp() << std::endl;
    }

    hpx::future<void> stage_future = receive_future.then([sim_unit](auto&& f) {
        f.get();  // check for exceptions

        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Starting work after receive" << std::endl;
        }

        sim_unit->discretization.mesh.CallForEachDistributedBoundary(
            [sim_unit](auto& dbound) { Problem::distributed_boundary_kernel(sim_unit->stepper, dbound); });

        sim_unit->discretization.mesh.CallForEachElement([sim_unit](auto& elt) {
            auto& state = elt.data.state[sim_unit->stepper.GetStage()];

            state.solution = elt.ApplyMinv(state.rhs);

            sim_unit->stepper.UpdateState(elt);
        });

        ++(sim_unit->stepper);

        if (sim_unit->writer.WritingVerboseLog()) {
            sim_unit->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
        }
    });

    if (SWE::PostProcessing::wetting_drying) {
        stage_future = stage_future.then([sim_unit](auto&& f) {
            f.get();  // check for exceptions

            sim_unit->discretization.mesh.CallForEachElement(
                [sim_unit](auto& elt) { Problem::wetting_drying_kernel(sim_unit->stepper, elt); });
        });
    }

    if (SWE::PostProcessing::slope_limiting) {
        stage_future = stage_future.then([sim_unit](auto&& f) {
            f.get();  // check for exceptions

            return CS_slope_limiter_hpx(sim_unit, CommTypes::baryctr_state);
        });
    }

    return stage_future.then([sim_unit](auto&& f) {
        f.get();  // check for exceptions

        sim_unit->discretization.mesh.CallForEachElement([sim_unit](auto& elt) {
            bool nan_found = SWE::scrutinize_solution(sim_unit->stepper, elt);

            if (nan_found)
                hpx::terminate();
        });
    });
}
}
}

#endif