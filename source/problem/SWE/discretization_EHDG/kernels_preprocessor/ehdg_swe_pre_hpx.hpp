#ifndef EHDG_SWE_PRE_HPX_HPP
#define EHDG_SWE_PRE_HPX_HPP

#include "problem/SWE/problem_preprocessor/swe_pre_init_data.hpp"

namespace SWE {
namespace EHDG {
template <typename HPXSimUnitType>
auto Problem::preprocessor_hpx(HPXSimUnitType* sim_unit) {
    initialize_data_parallel_pre_send(sim_unit->discretization.mesh, sim_unit->problem_input, CommTypes::baryctr_coord);

    hpx::future<void> receive_future =
        sim_unit->communicator.ReceiveAll(CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    sim_unit->communicator.SendAll(CommTypes::baryctr_coord, sim_unit->stepper.GetTimestamp());

    return receive_future.then([sim_unit](auto&&) {
        initialize_data_parallel_post_receive(sim_unit->discretization.mesh, CommTypes::baryctr_coord);

        Problem::initialize_global_problem(sim_unit->discretization);
    });
}
}
}

#endif