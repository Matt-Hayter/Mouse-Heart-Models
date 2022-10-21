#ifndef SIMULATION_CONFIGURATION_HPP
#define SIMULATION_CONFIGURATION_HPP

/*
Alter simulation_configuration to fit the required experiment.
*/

#include <string>

//Configuration class
struct simulation_configuration
{
    const std::string simulation_name = "RA"; //Name for simulation outputs

    const bool take_measurements = true; //Takes measurements of both the SAN and atrial action potentials
    const double dt = 0.001; //[ms] Time step for the simulation's solver
    const double output_dt = 0.01; //[ms] Intervals at which data is outputted. Must be divisible by dt
    
    const double temp = 310; //[K] Temperature for the simulation
    const int total_time = 3000; //[ms] Total time for the simulation

    //For SAN pacemaker - atrial myocyte junction:
    const bool homogeneous_coupling = false; //True: Homogeneous coupling conductance False: Individual coupling conductances
    double g_j_homogeneous = 1; //[nS] Homogeneous coupling conductance
    double g_j_atrial = 1.5; //[nS] Coupling for current flow in/out of atrial myocyte
    double g_j_SAN = 0; //[nS] Coupling for current flow in/out of SAN pacemaker
};

#endif