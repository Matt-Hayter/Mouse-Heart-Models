#ifndef SAN_CELL_HPP
#define SAN_CELL_HPP

/*
SAN pacemaking cell class
C++ implemented, coupled version of Wei's 2012 model
*/

#include "base_cell.hpp"
#include "simulation_configuration.hpp"
#include "SAN_cell_configuration.hpp"

#include <iostream>
#include <fstream>

//Configuration class
class SAN_cell : public base_cell
{
    protected:
        struct state_variables
        {
            double v;                       //1       
            double dst;                     //2       
            double fst;                     //3       
            double dt ;                     //4       
            double ft;                      //5       
            double ikr_act;                 //6           
            double ikr_act_f;               //7               
            double ikr_act_s;               //8               
            double ikr_inact;               //9               
            double iks_act;                 //10           
            double fl12;                    //11           
            double dl12;                    //12           
            double fl13;                    //13           
            double dl13;                    //14           
            double r;                       //15       
            double m_ttxr;                  //16           
            double h_ttxr;                  //17           
            double j_ttxr;                  //18           
            double m_ttxs;                  //19           
            double h_ttxs;                  //20           
            double j_ttxs;                  //21           
            double y_1;                     //22       
            double carel;                   //23           
            double caup;                    //24           
            double casub;                   //25           
            double Ftc;                     //26       
            double Ftmc;                    //27           
            double Ftmm;                    //28           
            double Fcms;                    //29           
            double Fcmi;                    //30           
            double Fcq;                     //31       
            double cai;                     //32       
            double q;                       //33       
            double fca;                     //34       
            double nai;                     //35       
            double ki;                      //36       
            double resting;                 //37           
            double open;                    //38           
            double inactivated;             //39               
            double resting_inactivated;     //40                       
            double y_2;                     //41       
            double y_4;                     //42    

            const int number_of_states{42}; //Stores the total number of state variables
        };
        struct outputs
        {
            double ist;                     //1 [pA] for all
            double ib;                      //2
            double ik1;                     //3
            double icat;                    //4
            double ikr;                     //5
            double iks;                     //6
            double ical12;                  //7
            double ical13;                  //8
            double ina_ttxs;                //9
            double ina_ttxr;                //10
            double ih;                      //11
            double ito;                     //12
            double isus;                    //13
            double inak;                    //14
            double inaca;                   //15
            double ina;                     //16
            double I_inwardcurrent;         //17
            double Itot;                    //18
            double I_j;                     //19 Junction current
        };

        //Instantiate states and configs
        const simulation_configuration sim_config;
        const SAN_configuration SAN_config;
        state_variables states; //States in the state vector
        state_variables states_dot; //Rates of states int the state vector
        outputs dependents;

        //Constants
        const double R = 8.31447;
        const double F = 96.4845;
        //Cell parameters
        const double capacitance = 0.025; //[nF]
        const double vcell = 3.0;
        const double l_cell = 66.3767257;
        const double r_cell = 3.792956;
        const double vrel = 0.0036;
        const double vsub = 0.03328117;
        const double vup = 0.0348;
        const double vi = 1.34671883;
        const double Mgi = 2.5;
        const double nao = 140.0;
        const double cao = 1.8;
        const double ko = 5.4;
        const double gst = 0.00006;
        const double eist = 17.0;
        const double gbna = 0.0001215;
        const double gbca = 0.000015;
        const double gbk = 0.0000025;
        const double gk1 = 0.229*0.0039228*0.9;
        const double gks = 0.000299;
        const double ecal = 47.0;
        const double kmfca = 0.00035;
        const double alpha_fca = 0.021;
        const double ecat = 45.0;
        const double enattxr = 41.5761;
        const double gsus = 0.00039060;
        const double inakmax = 0.7389*1.85*0.077;
        const double kmnap = 14.0;
        const double kmkp = 1.4;
        const double K1ni = 395.3;
        const double K1no = 1628;
        const double K2ni = 2.289;
        const double K2no = 561.4;
        const double K3ni = 26.44;
        const double K3no = 4.663;
        const double Kci = 0.0207;
        const double Kco = 3.663;
        const double Kcni = 26.44;
        const double Qci = 0.1369;
        const double Qco = 0.0;
        const double Qn = 0.4315;
        const double tdifca = 0.04;
        const double Prel = 2.5;
        const double Krel = 0.0015;
        const double nrel = 2.0;
        const double Kup = 0.0006;
        const double nup = 1.0;
        const double Ttr = 40.0;
        const double ConcTC = 0.031;
        const double ConcTMC = 0.062;
        const double kfTC = 88.8;
        const double kfTMC = 237.7;
        const double kbTC = 0.446;
        const double kbTMC = 0.00751;
        const double kfTMM = 2.277;
        const double kbTMM = 0.751;
        const double ConcCM = 0.045;
        const double kfCM = 237.7;
        const double kbCM = 0.542;
        const double ConcCQ = 10.0;
        const double kfCQ = 0.534;
        const double kbCQ = 0.445;
        const double koca = 10.0;
        const double kom = 0.06;
        const double kica = 0.5;
        const double kim = 0.005;
        const double eca50sr = 0.45;
        const double maxsr = 15.0;
        const double minsr = 1.0;
        const double hsrr = 2.5;
        const double pumphill = 2.0;
        const double gna_ttxs = 0.1*5.925e-05;
        const double gna_ttxr = 0.1*5.925e-05;
        //gCal*0.8 by wei
        const double gcal12 = 0.5*0.0010*4.0*1.5;
        const double gcal13 = 0.5*0.0030*4.0*1.5;
        //change by wei:
        const double gcat = 1.1*0.75*0.01862;
        const double ghk1 = 0.0015;
        const double ghk2 = 0.000728;
        const double ghk4 = 0.0025;
        const double ghna1 = 0.001;
        const double ghna2 = 0.000858;
        const double ghna4 = 0.00363;
        const double gkr = 0.8*0.002955;
        const double gto = 0.00492;
        const double kNaCa = 5.5;
        const double Pup = 0.04;
        const double ks = 1300000.0;
        const double Pup2 = 0.04;
        const double pumpkmf = 0.00008;
        const double pumpkmr = 4.5;

        double Vold{}; //Previous time step
        double dvdtold{};
        double top_slope_prev{}; //Previous beat
        double apd_start_prev{};

        int start_output{};
        double cycle_start;
        
        double min_potential;
        double tmin_potential;
        double max_potential;
        double dvdtmax;
        double apd_start;
        double ddr;
        double top;
        double top_slope;
        double apd30;
        double apd50;
        double apd90;
        double cycle_length;

        FILE *output_currents_file;
        FILE *output_measurements_file;

        double I_j; //[nA] Junction current

    public: 
        SAN_cell() {}
        ~SAN_cell() {}
        const double &get_potential() {return states.v;}
        void set_I_j(const double &junc_current) {I_j = junc_current*1e-3;} //[pA->nA]
        void output_config(std::ofstream &);
        void create_data_file();
        void create_measurement_file();
        void set_initial_states();
        void store_variables();
        void ODEs(const double &);
        void Euler_method();
        void output_data(const double &);
        void measurements(const double &, int &);
        void output_measurements(const int &);
        void output_final_states();
};

#endif