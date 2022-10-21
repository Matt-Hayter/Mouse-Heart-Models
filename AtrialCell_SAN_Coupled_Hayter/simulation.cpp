#include "simulation.hpp"

//Calculates the junction current between the SAN pacemaker and atrial myocyte
void simulation::process_junction_current(std::vector<base_cell*> &cells)
{
    double g_j_array[2]; //SAN, atrial
    double delta_V = cells[1]->get_potential() - cells[0]->get_potential(); //Atrial - SAN
    //Current = conductance * (V_atrial - V_SAN)
    if(sim_config.homogeneous_coupling == false) { //Isotropic currents
        g_j_array[0] = sim_config.g_j_SAN * delta_V; //[pA] SAN
        g_j_array[1] = sim_config.g_j_atrial * delta_V; //[pA] atrial
    } else {
        g_j_array[0,1] = sim_config.g_j_homogeneous * delta_V; //[pA]
    }
    for(int i{}; i < 2; i++) { //Set currents
        cells[i]->set_I_j(g_j_array[i]);
    }
}

const unsigned long long simulation::find_total_steps()
{
    const unsigned long long total_steps = sim_config.total_time/sim_config.dt; //Max number of steps for the simulation
    return total_steps;
}

//Prepares the simulation prior to solving the ODEs
void simulation::setup_simulation(std::vector<base_cell*> &cells)
{
    //Exception handling, ensuring valid user configuration
    try{
        output_total_config(cells);
        validate_sim_config();
        //Setup cells
        for(auto cell : cells) {
            cell->create_data_file();
            cell->create_measurement_file();
            cell->set_initial_states();
        }
        std::cout << "Beginning " << sim_config.simulation_name << " simulation..." << '\n' << std::endl;
    } catch(const std::exception &error) {
        std::cout << error.what() << std::endl; //Print error message associated with exception
        exit(EXIT_FAILURE);
    }
}

//Outputs simulation's configuration.
void simulation::output_total_config(const std::vector<base_cell*> &cells)
{
    //Simulation info
    std::ofstream config_file{".\\output_folder\\" + sim_config.simulation_name + "_config.txt"};
    if(config_file.fail()) throw std::runtime_error("output_folder not found in this directory"); //Check folder is present
    config_file << "Experiment type: " << "SAN coupled to atrial cell" << '\n' << "\n";
    config_file << "Simulation name: " << sim_config.simulation_name << '\n';
    config_file << "Solver time step(ms): " << sim_config.dt << '\n';
    config_file << "Output time step(ms): " << sim_config.output_dt << '\n';
    config_file << "Temperature(K): " << sim_config.temp << '\n';
    if(sim_config.homogeneous_coupling == false) {
        config_file << "Inomogeneous coupling: " << sim_config.g_j_SAN << '\n';
        config_file << "\tSAN g_j: " << sim_config.g_j_SAN << '\n';
        config_file << "\tAtrial g_j: " << sim_config.g_j_atrial << '\n' << '\n';
    } else {
        config_file << "Homogeneous coupling: " << sim_config.g_j_SAN << '\n';
        config_file << "\tg_j: " << sim_config.g_j_homogeneous << '\n' << '\n';
    }

    //Ouput configs for both cells
    for (auto cell : cells) {
        cell->output_config(config_file);
    }
}

//Validates the user's core configuration setup
void simulation::validate_sim_config()
{
    //Check time-steps
    if(fabs(remainder(sim_config.output_dt, sim_config.dt)) > 1e-12) 
        throw std::logic_error("Invalid 'dt' and 'output_dt' configured. Check output_dt is divisible by dt.");
    if(sim_config.output_dt <= 0 || sim_config.dt <= 0) throw std::logic_error("Invalid 'dt' and 'output_dt' configured. Check that output_dt and dt are greater than 0.");
}

//Update simulation time
void simulation::update_time(const int &step)
{
    time = step * sim_config.dt; //Update time
}

void simulation::end_simulation(std::vector<base_cell*> &cells)
{
    try{
        //Ouput final states for each cell
        for(auto cell : cells) {
            cell->output_final_states();
        }
    } catch(const std::exception &error) {
        std::cout << error.what() << std::endl; //Print error message associated with exception
    }
    std::cout << '\n' << "Simulation Completed";
}