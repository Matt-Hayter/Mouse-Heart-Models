#include "SAN_cell.hpp"
#include <cmath>

void SAN_cell::output_config(std::ofstream &config_file)
{
    //SAN cell configuration:
    config_file << "-----------------------SAN Configuration-------------------------------" << '\n' << '\n';
    config_file << "Initial states file name: " << SAN_config.initial_filename << '\n';
    config_file << "Final states file name: " << sim_config.simulation_name << "_SAN_final_states.dat" << '\n' << "\n";
    config_file << "INa Block: " << SAN_config.INa_Block << '\n' << "\n";
}

void SAN_cell::create_data_file()
{
    //File to output currents, in a formatted string
    std::string file_path{".\\output_folder\\" + sim_config.simulation_name + "_SAN_data.dat"};
    output_currents_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    //Output and format column labels
    fprintf(output_currents_file, "\t%-17s%-25s%-15s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s\n",
                "Time(ms)", "Membrane Potential(mV)", "dV/dt(V/s)", "I_Na", "INa_ttxr", "INa_ttxs", "ICaT", "ICaL12", "ICaL13",
                "IKr", "IKs", "Ito", "INaCa", "Ih", "INaK", "Ist", "IK1", "Isus", "Ib", "I_inward", "I_tot"
    );
}

void SAN_cell::create_measurement_file()
{
    std::string file_path{".\\output_folder\\" + sim_config.simulation_name + "_SAN_measurements.dat"};
    output_measurements_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    fprintf(output_measurements_file, "\t%-15s%-20s%-20s%-15s%-15s%-20s%-20s%-20s%-15s%-17s%-15s\n", "Beat Number", "AP Start Time(ms)", "Max Potential(mV)", "MDP(mV)", "Max dV/dt(V/s)", "APD_30(ms)", "APD_50(ms)", "APD_90(ms)", "APD_50/APD_90", "Cycle Length(ms)", "Cycle Frequency(Hz)");
}

//Parse initial states datafile, validate and add to state vector
void SAN_cell::set_initial_states()
{
    double* state_ptr = &states.v; //Pointer to first state

    std::ifstream initial_states_file;
    initial_states_file.open("./input_folder/SAN_inputs/" + SAN_config.initial_filename);
    if(initial_states_file.fail()) throw std::runtime_error("SAN initial states file not found, under ./input_folder/SAN_inputs/...");
    std::string line;
    //Loop for reading in states
    for(int i{}; i < states.number_of_states; i++) {
        if(getline(initial_states_file, line).fail() || (i + 1 == states.number_of_states && initial_states_file.peek() != -1)) { //Read in state and check number
            throw std::logic_error("Invalid SAN initial states. " + std::to_string(i) + " states detected in initial states file but expected "
                + std::to_string(states.number_of_states) + ". Check initial states file or number_of_states variable definition.");
        }
        *state_ptr = atof(line.c_str()); //Convert string to double
        *(++state_ptr); //Move to next state
    }
}

void SAN_cell::output_data(const double &time)
{
    if(fabs(remainder(time, sim_config.output_dt)) > 1e-12) return; //Ensuring data is only outputted on the configured intervals
    fprintf(output_currents_file, "\t%-17.10g%-25.10g%-15g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g\n",
                time, states.v, states_dot.v, dependents.ina, dependents.ina_ttxr, dependents.ina_ttxs, dependents.icat, dependents.ical12, 
                dependents.ical13, dependents.ikr, dependents.iks, dependents.ito, dependents.inaca, dependents.ih, dependents.inak,
                dependents.ist, dependents.ik1, dependents.isus, dependents.ib, dependents.I_inwardcurrent, dependents.Itot
    );
}

void SAN_cell::Euler_method()
{
    states.v += sim_config.dt*states_dot.v;
    states.dst += sim_config.dt*states_dot.dst;
    states.fst += sim_config.dt*states_dot.fst;
    states.dt += sim_config.dt*states_dot.dt ;
    states.ft += sim_config.dt*states_dot.ft;
    //states.ikr_act += sim_config.dt*states_dot.ikr_act;
    states.ikr_act_f += sim_config.dt*states_dot.ikr_act_f;
    states.ikr_act_s += sim_config.dt*states_dot.ikr_act_s;
    states.ikr_inact += sim_config.dt*states_dot.ikr_inact;
    states.iks_act += sim_config.dt*states_dot.iks_act;
    states.fl12 += sim_config.dt*states_dot.fl12;
    states.dl12 += sim_config.dt*states_dot.dl12;
    states.fl13 += sim_config.dt*states_dot.fl13;
    states.dl13 += sim_config.dt*states_dot.dl13;
    states.r += sim_config.dt*states_dot.r;
    states.m_ttxr += sim_config.dt*states_dot.m_ttxr;
    states.h_ttxr += sim_config.dt*states_dot.h_ttxr;
    states.j_ttxr += sim_config.dt*states_dot.j_ttxr;
    states.m_ttxs += sim_config.dt*states_dot.m_ttxs;
    states.h_ttxs += sim_config.dt*states_dot.h_ttxs;
    states.j_ttxs += sim_config.dt*states_dot.j_ttxs;
    states.y_1 += sim_config.dt*states_dot.y_1;
    states.carel += sim_config.dt*states_dot.carel;
    states.caup += sim_config.dt*states_dot.caup;
    states.casub += sim_config.dt*states_dot.casub;
    states.Ftc += sim_config.dt*states_dot.Ftc;
    states.Ftmc += sim_config.dt*states_dot.Ftmc;
    states.Ftmm += sim_config.dt*states_dot.Ftmm;
    states.Fcms += sim_config.dt*states_dot.Fcms;
    states.Fcmi += sim_config.dt*states_dot.Fcmi;
    states.Fcq += sim_config.dt*states_dot.Fcq;
    states.cai += sim_config.dt*states_dot.cai;
    states.q += sim_config.dt*states_dot.q;
    states.fca += sim_config.dt*states_dot.fca;
    states.nai -= sim_config.dt*states_dot.nai;
    states.ki -= sim_config.dt*states_dot.ki;
    states.resting += sim_config.dt*states_dot.resting;
    states.open += sim_config.dt*states_dot.open;
    states.inactivated += sim_config.dt*states_dot.inactivated;
    states.resting_inactivated += sim_config.dt*states_dot.resting_inactivated;
    states.y_2 += sim_config.dt*states_dot.y_2;
	states.y_4 += sim_config.dt*states_dot.y_4;
}

void SAN_cell::measurements(const double &time, int &current_beat)
{   
    double *dvdt = &states_dot.v; //For conciseness

    if(*dvdt>=0.0&&dvdtold<0.0) {
        //Now progress to next beat:
        current_beat += 1;
        std::cout << '\n' << "Generating AP " << current_beat;
        if(sim_config.take_measurements == true) { //Take measurements for previous beat
            //Finalise measurements for beat
            cycle_length = time - cycle_start;
            //Output data after full beat:
            if(current_beat > 1) { //Don't output measurements on first beat
                std::cout <<": " << std::endl; //For output formatting
                output_measurements(current_beat);
            }
            //Set params for AP:
            cycle_start = time; 
            start_output = 1;
            min_potential = Vold;
            tmin_potential = time;
            dvdtmax = -10000;
        } else {
            std::cout << std::endl; //For output formatting
        }
    }
    if(sim_config.take_measurements == true) { //If true, take measurements
        if(*dvdt>dvdtmax&&start_output>0) {
            dvdtmax = *dvdt;
            apd_start = time; //Start APD measurements at max dV/dt
        }

        if(dvdtold>0.0&&*dvdt<=0.0) {
            max_potential = Vold;
            top_slope = (max_potential-min_potential)/(time - tmin_potential);
        }

        if((dvdtold<=top_slope_prev)&&(*dvdt>top_slope_prev)) {
            top = Vold;
            ddr = (Vold - min_potential)/(time - tmin_potential);
        }
        if(states.v<=0.7*max_potential+0.3*min_potential&&Vold>0.7*max_potential+0.3*min_potential) {
            apd30 = time - apd_start;
        }
        if(states.v<=0.5*max_potential+0.5*min_potential&&Vold>0.5*max_potential+0.5*min_potential) {
            apd50 = time - apd_start;
        }

        if(states.v<=0.1*max_potential+0.9*min_potential&&Vold>0.1*max_potential+0.9*min_potential)
        {
            if(apd_start>0.0) {
                apd90 = time - apd_start;
                //Store values for comparison on next beat
                top_slope_prev = top_slope;
            }
        }
    }
}

void SAN_cell::output_measurements(const int &current_beat)
{
    double norm_apd50 = apd50/apd90;
    double cycle_freq = 1/(cycle_length/1000);
    fprintf(output_measurements_file, "\t%-15i%-20.10g%-20.6g%-15.6g%-15.5g%-20.6g%-20.6g%-20.6g%-15.6g%-17.6g%-15.6g\n", current_beat, cycle_start, max_potential, min_potential, dvdtmax, apd30, apd50, apd90, norm_apd50, cycle_length, cycle_freq);
    std::cout << "\tSAN cell measurements:" << '\n';
    std::cout << "\t\tSAN Cycle start time(ms): " << cycle_start << '\n';
    std::cout << "\t\tSAN Max potential(mV): " << max_potential << '\n';
    std::cout << "\t\tSAN MDP(mV): " << min_potential << '\n';
    std::cout << "\t\tSAN Max dV/dt(V/s): " << dvdtmax << '\n';
    std::cout << "\t\tSAN APD30(ms): " << apd30 << '\n';
    std::cout << "\t\tSAN APD50(ms): " << apd50 << '\n';
    std::cout << "\t\tSAN APD90(ms): " << apd90 << '\n';
    std::cout << "\t\tSAN Cycle length(ms): " << cycle_length << '\n';
    std::cout << "\t\tSAN Cycle frequency(Hz): " << cycle_freq << '\n';
    std::cout << "\t\tSAN DDR: " << ddr << '\n';
    std::cout << "\t\tSAN TOP: " << top << std::endl;
}

void SAN_cell::store_variables()
{
    Vold = states.v;
    dvdtold = states_dot.v;
}

//Output the final state variables
void SAN_cell::output_final_states()
{
    const double *state_ptr = &states.v; //Pointer to the first state

    std::ofstream final_states_file{".\\output_folder\\" + sim_config.simulation_name + "_SAN_final_states.dat"};
    if(final_states_file.fail()) throw std::runtime_error{"output_folder was not found at the end of the simulation."};
    final_states_file.precision(20);
    for(int i{}; i < states.number_of_states; i++) { //Output each state
        final_states_file << *state_ptr << std::endl;
        *(++state_ptr);
    }
}