/*
Matt Hayter

Simplified atrial myocyte model coupled to a SAN pacemaker cell model.
*/

#include "simulation_configuration.hpp"
#include "simulation.hpp"

#include "SAN_cell_configuration.hpp"
#include "SAN_cell.hpp"

#include "atrial_cell_configuration.hpp"
#include "atrial_cell.hpp"

int main()
{
    simulation sim; //Instantiate simulation and cells
    std::vector<base_cell*> cells; //Create cell vector and append desired cells
    cells.push_back(new SAN_cell);
    cells.push_back(new atrial_cell); //Add additional cells here

    sim.setup_simulation(cells);

    const int total_steps = sim.find_total_steps(); //Number of steps for the simulation

    for(int step{1}; step <= total_steps; step++) {
        sim.process_junction_current(cells);
        for(auto cell : cells) {
            cell->store_variables(); //Store variables for next step
            cell->ODEs(sim.time);
            cell->Euler_method(); //Progress state variables to next step
        }
        sim.update_time(step); //Progress simulation time step
        for(auto cell : cells) {
            cell->measurements(sim.time, sim.current_beat);
            cell->output_data(sim.time);
        }
    }
    sim.end_simulation(cells);
    
    delete cells[0], cells[1];
    return 0;
}