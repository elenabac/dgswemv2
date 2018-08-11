#ifndef RKDG_SWE_PRE_SERIAL_HPP
#define RKDG_SWE_PRE_SERIAL_HPP

#include "rkdg_swe_pre_init_data.hpp"

namespace SWE {
namespace RKDG {
void Problem::serial_preprocessor_kernel(ProblemDiscretizationType& discretization,
                                         const ProblemInputType& problem_specific_input) {
    Problem::initialize_data_kernel(discretization.mesh, problem_specific_input);
}
}
}

#endif