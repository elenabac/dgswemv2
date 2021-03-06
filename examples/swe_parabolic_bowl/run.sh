#!/bin/bash

DGSWEMV2_ROOT=/home/bachini/C++/dgswemv2
DGSWEMV2_BUILD_DIR=$DGSWEMV2_ROOT/build

# This is a walk through of the dgswemv2 workflow.  This script
# recreates the parabolic bowl solution in:
# Bunya et al., "A Wetting
# and Drying Treatment for the Runge-Kutta Discontinuous Galerkin
# Solution to the Shallow Water Equations", CMAME (2009).

# Throughout the simulation, we will rely on the following environment
# variables.  DGSWEMV2_ROOT should be set to the dgswemv2 repository
# path.
if [[ ! -d $DGSWEMV2_ROOT ]]; then
    echo "Error: Could not find DGSWEMV2_ROOT: ${DGSWEMV2_ROOT}"
    echo "       Please set the environment variable DGSWEMV2_ROOT to point at the"
    echo "       dgswemv2 repository."
    exit 1
fi

#We use internally use DGSWEMV2_ROOT_, which we set to an absolute path.
if [[ ! "${DGSWEMV2_ROOT}" = /* ]]; then
    DGSWEMV2_ROOT_="$(pwd)/${DGSWEMV2_ROOT}"
else
    DGSWEMV2_ROOT_=$DGSWEMV2_ROOT
fi

# BUILD_DIR should be set to build directory for dgswemv2
if [[ ! -d $DGSWEMV2_BUILD_DIR &&
          ! -d $DGSWEMV2_ROOT/$DGSWEMV2_BUILD_DIR ]]; then
    echo "Error: Could not find DGSWEMV2_BUILD_DIR: ${DGSWEMV2_BUILD_DIR}"
    echo "       Please set the environment variable DGSWEMV2_BUILD_DIR to point at the"
    echo "       build directory. Note that either absolute or relative paths will work"
    exit 1
fi

#Similarly, we internally use DGSWEMV2_BUILD_DIR_, which we set to an absolute path
if [[ "${DGSWEMV2_BUILD_DIR}" = /* ]]; then
    DGSWEMV2_BUILD_DIR_=$DGSWEMV2_BUILD_DIR
else
    if [[ -d $DGSWEMV2_BUILD_DIR ]]; then
        DGSWEMV2_BUILD_DIR_="$(pwd)/${DGSWEMV2_BUILD_DIR}"
    else
        DGSWEMV2_BUILD_DIR_="${DGSWEMV2_ROOT_}/${DGSWEMV2_BUILD_DIR}"
    fi
fi

PARABOLIC_BOWL_DIR=$DGSWEMV2_ROOT_/examples/swe_parabolic_bowl

###################################################################################
echo "Running the Simulation"
###################################################################################

#To run the simulation, similar to the mesh generation phase, we need
#to copy the initial conditions, and recompile the dgswemv2 target.
cd $DGSWEMV2_ROOT_/source/problem/SWE/problem_function_files
mv swe_initial_condition_functions.hpp swe_initial_condition_functions.hpp.tmp
cp $PARABOLIC_BOWL_DIR/parabolic_bowl_initial_condition_functions.hpp swe_initial_condition_functions.hpp
mv swe_true_solution_functions.hpp swe_true_solution_functions.hpp.tmp
cp $PARABOLIC_BOWL_DIR/parabolic_bowl_true_solution_functions.hpp swe_true_solution_functions.hpp

#compute L2Residual
MAIN_DIR="$DGSWEMV2_ROOT_/source"
cp ${MAIN_DIR}/dgswemv2-serial.cpp ${MAIN_DIR}/dgswemv2-serial.cpp.tmp
sed -i.tmp '/return 0/i\
simulation->ComputeL2Residual();' ${MAIN_DIR}/dgswemv2-serial.cpp

#Now we compile the dgswemv2-serial target
cd $DGSWEMV2_BUILD_DIR_
make dgswemv2-serial

#We are now ready to run the simulation.
#Begin by going to the parabolic bowl directory, and making the output directory
cd $PARABOLIC_BOWL_DIR/input_files

mkdir -p output

#the dgswemv2-serial target requires a single argument of the input file
#In general, paths can be relative or absolute. Here we use relative paths.
$DGSWEMV2_BUILD_DIR_/source/dgswemv2-serial dgswemv2_input.15

#These configuration will write vtu outputs to the output directory.

###################################################################################
echo "Cleaning up"
###################################################################################
cd $DGSWEMV2_ROOT_/source/problem/SWE/problem_function_files
mv swe_initial_condition_functions.hpp.tmp swe_initial_condition_functions.hpp
mv swe_true_solution_functions.hpp.tmp swe_true_solution_functions.hpp
mv ${MAIN_DIR}/dgswemv2-serial.cpp.tmp ${MAIN_DIR}/dgswemv2-serial.cpp

cd $DGSWEMV2_ROOT_/mesh_generators
if [[ -f bathymetry.hpp.tmp ]]; then
    mv bathymetry.hpp.tmp bathymetry.hpp
fi
