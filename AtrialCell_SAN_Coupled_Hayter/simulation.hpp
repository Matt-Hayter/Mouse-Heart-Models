#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "atrial_cell.hpp"
#include "SAN_cell.hpp"
#include "simulation_configuration.hpp"

//Class for simulation as a whole

struct simulation
{
    double time{}; //Current time for simulation
    int current_beat{}; //Tracks the number of the current SAN AP being generated
    bool AP_start; //True for the first time-step of outward SAN junction current for each AP
    bool take_atrial_measurements; //Required to instruct the atrial cell when to start/stop taking measurements

    //Declare simulation configuration
    simulation_configuration sim_config;

    const unsigned long long find_total_steps();
    void setup_simulation(std::vector<base_cell*> &);
    void output_total_config(const std::vector<base_cell*> &);
    void validate_sim_config();
    void output_sim_measurements();
    void update_time(const int &);
    void end_simulation(std::vector<base_cell*> &);
    void process_junction_current(std::vector<base_cell*> &);
};

#endif