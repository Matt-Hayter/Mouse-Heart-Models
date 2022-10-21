#ifndef SAN_CELL_CONFIGURATION_HPP
#define SAN_CELL_CONFIGURATION_HPP

/*
Alter SAN_cell_configuration to fit the required experiment.

All other SAN cell parameters can be configured from within SAN_cell's member data, or within related ODE functions
*/

#include <string>

//Configuration class
struct SAN_configuration
{
    const std::string initial_filename = "SAN_final_states.dat"; //File containing initial states to be parsed. Must be within "input_folder/SAN_inputs"

    //Channel blocks - Fraction of remaining current
    const double INa_Block = 1;
    const double ICaL_Block = 1;
};

#endif