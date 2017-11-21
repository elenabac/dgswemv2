#ifndef OMPI_SIMULATION_HPP
#define OMPI_SIMULATION_HPP

#include "../general_definitions.hpp"

#include "../preprocessor/input_parameters.hpp"
#include "../preprocessor/initialize_mesh.hpp"
#include "../communication/ompi_communicator.hpp"
#include "writer.hpp"

#include "utilities/file_exists.hpp"

template <typename ProblemType>
class OMPISimulationUnit {
  private:
    InputParameters<typename ProblemType::ProblemInputType> input;

    Stepper stepper;
    Writer<ProblemType> writer;
    typename ProblemType::ProblemMeshType mesh;
    OMPICommunicator communicator;

  public:
    OMPISimulationUnit()
        : input(),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          mesh(this->input.polynomial_order) {}
    OMPISimulationUnit(const std::string& input_string, const uint locality_id, const uint submesh_id)
        : input(input_string, locality_id, submesh_id),
          stepper(this->input.rk.nstages, this->input.rk.order, this->input.dt),
          writer(input, locality_id, submesh_id),
          mesh(this->input.polynomial_order),
          communicator(this->input.mesh_file_name.substr(0, this->input.mesh_file_name.find_last_of('.')) + ".dbmd",
                       locality_id,
                       submesh_id) {
        this->input.ReadMesh();

        this->mesh.SetMeshName(this->input.mesh_data.mesh_name);

        if (this->writer.WritingLog()) {
            this->writer.StartLog();

            this->writer.GetLogFile() << "Starting simulation with p=" << input.polynomial_order << " for "
                                      << input.mesh_data.mesh_name << " mesh" << std::endl << std::endl;
        }

        initialize_mesh<ProblemType, OMPICommunicator>(
            this->mesh, this->input.mesh_data, this->communicator, this->input.problem_input, this->writer);
    }

    void Launch();

    void ExchangeData();
    void PreReceiveStage();
    void PostReceiveStage();
    void WaitAllSends();

    void Step();

    void ResidualL2();
};

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Launch() {
    if (this->writer.WritingLog()) {
        this->writer.GetLogFile() << std::endl << "Launching Simulation!" << std::endl << std::endl;
    }

    if (this->writer.WritingOutput()) {
        this->writer.WriteFirstStep(this->stepper, this->mesh);
    }

    this->communicator.InitializeCommunication();

    uint n_stages = this->stepper.get_num_stages();

    auto resize_data_container = [n_stages](auto& elt) { elt.data.resize(n_stages); };

    this->mesh.CallForEachElement(resize_data_container);
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::ExchangeData() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Current (time, stage): (" << this->stepper.get_t_at_curr_stage() << ','
                                  << this->stepper.get_stage() << ')' << std::endl;
    }

    auto distributed_boundary_send_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_send_kernel(this->stepper, dbound);
    };

    this->communicator.ReceiveAll(this->stepper.get_timestamp());

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_send_kernel);

    this->communicator.SendAll(this->stepper.get_timestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PreReceiveStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work before receive" << std::endl;
    }

    auto volume_kernel = [this](auto& elt) { ProblemType::volume_kernel(this->stepper, elt); };

    auto source_kernel = [this](auto& elt) { ProblemType::source_kernel(this->stepper, elt); };

    auto interface_kernel = [this](auto& intface) { ProblemType::interface_kernel(this->stepper, intface); };

    auto boundary_kernel = [this](auto& bound) { ProblemType::boundary_kernel(this->stepper, bound); };

    this->mesh.CallForEachElement(volume_kernel);

    this->mesh.CallForEachElement(source_kernel);

    this->mesh.CallForEachInterface(interface_kernel);

    this->mesh.CallForEachBoundary(boundary_kernel);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work before receive" << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::PostReceiveStage() {
    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting to wait on receive with timestamp: " << this->stepper.get_timestamp()
                                  << std::endl;
    }

    this->communicator.WaitAllReceives(this->stepper.get_timestamp());

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Starting work after receive" << std::endl;
    }

    auto distributed_boundary_kernel = [this](auto& dbound) {
        ProblemType::distributed_boundary_kernel(this->stepper, dbound);
    };

    auto update_kernel = [this](auto& elt) { ProblemType::update_kernel(this->stepper, elt); };

    auto scrutinize_solution_kernel = [this](auto& elt) {
        ProblemType::scrutinize_solution_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachDistributedBoundary(distributed_boundary_kernel);

    this->mesh.CallForEachElement(update_kernel);

    this->mesh.CallForEachElement(scrutinize_solution_kernel);

    ++(this->stepper);

    if (this->writer.WritingVerboseLog()) {
        this->writer.GetLogFile() << "Finished work after receive" << std::endl << std::endl;
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::WaitAllSends() {
    this->communicator.WaitAllSends(this->stepper.get_timestamp());
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::Step() {
    auto swap_states_kernel = [this](auto& elt) { ProblemType::swap_states_kernel(this->stepper, elt); };

    this->mesh.CallForEachElement(swap_states_kernel);

    if (this->writer.WritingOutput()) {
        this->writer.WriteOutput(this->stepper, mesh);
    }
}

template <typename ProblemType>
void OMPISimulationUnit<ProblemType>::ResidualL2() {
    double residual_L2 = 0;

    auto compute_residual_L2_kernel = [this, &residual_L2](auto& elt) {
        residual_L2 += ProblemType::compute_residual_L2_kernel(this->stepper, elt);
    };

    this->mesh.CallForEachElement(compute_residual_L2_kernel);

    this->writer.GetLogFile() << "residual inner product: " << residual_L2 << std::endl;
}

template <typename ProblemType>
class OMPISimulation {
  private:
    uint n_steps;
    uint n_stages;

    std::vector<std::unique_ptr<OMPISimulationUnit<ProblemType>>> simulation_units;

  public:
    OMPISimulation() = default;
    OMPISimulation(const std::string& input_string) {
        int locality_id;
        MPI_Comm_rank(MPI_COMM_WORLD, &locality_id);

        InputParameters<typename ProblemType::ProblemInputType> input(input_string);

        this->n_steps = (uint)std::ceil(input.T_end / input.dt);
        this->n_stages = input.rk.nstages;

        std::string submesh_file_prefix = input.mesh_file_name.substr(0, input.mesh_file_name.find_last_of('.')) + "_" +
                                          std::to_string(locality_id) + '_';
        std::string submesh_file_postfix =
            input.mesh_file_name.substr(input.mesh_file_name.find_last_of('.'), input.mesh_file_name.size());

        uint submesh_id = 0;

        while (Utilities::file_exists(submesh_file_prefix + std::to_string(submesh_id) + submesh_file_postfix)) {
            this->simulation_units.emplace_back(
                new OMPISimulationUnit<ProblemType>(input_string, locality_id, submesh_id));

            ++submesh_id;
        }
    }

    void Run();
};

template <typename ProblemType>
void OMPISimulation<ProblemType>::Run() {
#pragma omp parallel
    {
        uint n_threads, thread_id, sim_per_thread, begin_sim_id, end_sim_id;

        n_threads = (uint)omp_get_num_threads();
        thread_id = (uint)omp_get_thread_num();

        sim_per_thread = (this->simulation_units.size() + n_threads - 1) / n_threads;

        begin_sim_id = sim_per_thread * thread_id;
        end_sim_id = std::min(sim_per_thread * (thread_id + 1), (uint) this->simulation_units.size());

        for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
            this->simulation_units[sim_unit_id]->Launch();
        }

        for (uint step = 1; step <= this->n_steps; step++) {
            for (uint stage = 0; stage < this->n_stages; stage++) {
                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->ExchangeData();
                }

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->PreReceiveStage();
                }

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->PostReceiveStage();
                }

                for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                    this->simulation_units[sim_unit_id]->WaitAllSends();
                }
            }

            for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
                this->simulation_units[sim_unit_id]->Step();
            }
        }
#ifdef RESL2
        for (uint sim_unit_id = begin_sim_id; sim_unit_id < end_sim_id; sim_unit_id++) {
            this->simulation_units[sim_unit_id]->ResidualL2();
        }
#endif
    }
}

#endif
