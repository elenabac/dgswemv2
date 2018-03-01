#ifndef SWE_INPUTS_HPP
#define SWE_INPUTS_HPP

#include "../../../general_definitions.hpp"
#include <yaml-cpp/yaml.h>

namespace SWE {
enum class SphericalProjectionType {
    None,
    Enable
};

struct SphericalProjection {
    SphericalProjectionType type = SphericalProjectionType::None;
    double polar_o = 0.0;
    double azimuthal_o = 0.0;
    double R = 1000.0;
};

enum class InitialConditionsType {
    Constant,
    Function
};

struct InitialConditions {
    InitialConditionsType type = InitialConditionsType::Constant;
    double ze_initial = 0.;
    double qx_initial = 0.;
    double qy_initial = 0.;
};

enum class FunctionSourceType {
    None,
    Enable
};

struct FunctionSource {
    FunctionSourceType type = FunctionSourceType::None;
};

enum class BottomFrictionType {
    None,
    Chezy,
    Manning
};

struct BottomFriction {
    BottomFrictionType type = BottomFrictionType::None;
    double coefficient = 0.0;
    std::string manning_data_file;
};

enum class MeteoForcingType {
    None,
    Test
};

struct MeteoForcing {
    MeteoForcingType type = MeteoForcingType::None;
    std::string meteo_data_file;
    double frequency;
};

enum class TidalPotentialType {
    None,
    Test
};

struct TidalPotential {
    TidalPotentialType type = TidalPotentialType::None;
};

enum class CoriolisType {
    None,
    Enable
};

struct Coriolis {
    CoriolisType type = CoriolisType::None;
};

struct Inputs {
    Inputs() = default;
    Inputs(YAML::Node& swe_node);

    YAML::Node as_yaml_node();

    double g = 9.81;
    double h_o = 0.01;

    bool parse_input = false;

    SphericalProjection spherical_projection;

    InitialConditions initial_conditions;

    FunctionSource function_source;
    BottomFriction bottom_friction;
    MeteoForcing meteo_forcing;
    TidalPotential tidal_potential;
    Coriolis coriolis;
};
}

#endif