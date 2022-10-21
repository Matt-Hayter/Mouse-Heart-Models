/*
Matt Hayter

Simpified version of the H. Zhang et al. (2020) Mouse Atrial Cell Model, with ODEs 
solved using the implicit Euler method.

This version encorporates steady pacing of the model only.

02/2022
*/

#include "configuration.hpp"
#include "atrial_cell.hpp"
#include <chrono>

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    atrial_cell atrial_myocyte;
    const atrial_configuration &config = atrial_myocyte.get_config(); //Fetch config for use in main

    setup_simulation(atrial_myocyte, config);

    const int max_steps = (config.beats*config.BCL)/config.dt; //Max number of steps for the simulation

    for(int step{1}; step <= max_steps; step++) {
        atrial_myocyte.store_variables();
        atrial_myocyte.process_stimulus(atrial_myocyte.get_time()); //Generates stimulus
        
        atrial_myocyte.IP3_ODEs();
        atrial_myocyte.cam_ODEs();
        atrial_myocyte.camkii_ODEs();
        atrial_myocyte.bars_ODEs();
        atrial_myocyte.ecc_ODEs();

        //Progress to the next time step
        atrial_myocyte.Euler_method();
        atrial_myocyte.update_time(step);

        atrial_myocyte.measurements();

        atrial_myocyte.output_data();
        atrial_myocyte.output_measurements();
    }

    atrial_myocyte.end_simulation();

    auto t2 = std::chrono::high_resolution_clock::now();

    auto second_int = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    std::cout << '\n' << second_int.count();
    return 0;
}