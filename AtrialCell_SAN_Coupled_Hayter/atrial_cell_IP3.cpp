#include "atrial_cell.hpp"

void atrial_cell::IP3_ODEs(const double &time)
{
    //For conciseness
    const ip3s *IP3 = &states.ip3;
    eccs *ec = &states.ec;
    ip3s *IP3R = &states_dot.ip3;
    eccs *ecR = &states_dot.ec;

    double Ca_i = states.ec.Ca_i * 1e3;    //Cytosolic Ca concentration, [mM ->uM]

    if (time > atrial_config.IP3_TIME) {
        if (!atrial_config.IP3_dynamic_on) {
            ec->IP3_state = atrial_config.IP3_fixed; //No IP3 dynamic, set it to a fixed value
            ecR->IP3_state = 0; //So the state remains constant
        } else { //IP3 is dynamic, so give ODEs
            double L = atrial_config.ET_1; //Ligand stimulation concentration, [uM]
            double PIP2 = 4000; //PIP2 density, [um^2-1]
            double kf1 = 0.0003; //R1 forward rate constant, [uM-1 s-1]
            double Kd1 = 3e-5; //R1 dissociation constant, [uM]
            double kf2 = 2.75e-4; //[um^2/s]
            double Kd2 = 27500; //[um^2-1]
            double kf3 = 1; //[um^2/s]
            double kr3 = 0.001; //[s-1]
            double kf4 = 0.3; //[uM-1 s-1]
            double Kd4 = 3e-5; //[uM]
            double kf5 = 0.0004; //[s-1]
            double kf6 = 1; //[s-1]
            double kf7 = 0.15; //[s-1]
            double kf8 = 0.0167; //[uM-1 s-1]
            double kr8 = 0.0167; //[s-1]
            double kf9 = 0.0042; //[um^2/s]
            double kr9 = 1; //[s-1]
            double kf10 = 0.042; //[um^2/s]
            double kr10 = 1; //[s-1]
            double kf11 = 0.0334; //[uM-1 s-1]
            double Kd11 = 0.1; //[uM]
            double kf12 = 6; //[s-1]
            double kf13 = 6; //[s-1]
            double kf14 = 0.444; //[s-1]
            double Km14 = 19.8; //[uM]  
            double kf15 = 3.8; //[s-1]
            double Km15 = 5; //[uM]
            double kf16 = 1.25; //[s-1]
            double Vc = 2550; //[um^3]
            double Rpc = 4.61; //[um-1]

            double kr1 = kf1 * Kd1; //[s-1]
            double kr2 = kf2 * Kd2; //[s-1]
            double kr4 = kf4 * Kd4; //[s-1]
            double kr11 = kf11 * Kd11; //[s-1]
            double Cc = 1.00000 / (Vc * 602.200); //[uM]
            double Cp = 1.00000 / (Vc * Rpc); //[um^2-1]
            double Cpc = Cc / Cp; //[uM um^2]

            double J1 = 1e-3 * (kf1 * IP3->r * L - kr1 * IP3->Rl); //[um^2-1 ms-1]
            double J2 = 1e-3 * (kf2 * IP3->r * IP3->Gd - kr2 * IP3->Rg); //[um^2-1 ms-1]
            double J3 = 1e-3 * (kf3 * IP3->Rl * IP3->Gd - kr3 * IP3->Rlg); //[um^2-1 ms-1]
            double J4 = 1e-3 * (kf4 * L * IP3->Rg - kr4 * IP3->Rlg); //[um^2-1 ms-1]
            double J5 = 1e-3 * kf5 * IP3->Rlg; //[um^2-1 ms-1]
            double J6 = 1e-3 * kf6 * IP3->Rlg; //[um^2-1 ms-1]
            double J7 = 1e-3 * kf7 * IP3->Gt; //[um^2-1 ms-1]
            double J8 = 1e-3 * (kf8 * IP3->P * Ca_i - kr8 * IP3->Pc); //[um^2-1 ms-1]
            double J9 = 1e-3 * (kf9 * IP3->P * IP3->Gt - kr9 * IP3->Pg); //[um^2-1 ms-1]
            double J10 = 1e-3 * (kf10 * IP3->Pc * IP3->Gt - kr10 * IP3->Pcg); //[um^2-1 ms-1]
            double J11 = 1e-3 * (kf11 * IP3->Pg * Ca_i - kr11 * IP3->Pcg); //[um^2-1 ms-1]
            double J12 = 1e-3 * kf12 * IP3->Pcg; //[um^2-1 ms-1]
            double J13 = 1e-3 * kf13 * IP3->Pg; //[um^-2 * s^-1] -> [um^-2 * ms^-1]
            double J14 = 1e-3 * (kf14 * IP3->Pc * PIP2) / (Km14 / Cpc + PIP2); //[um^2-1 ms-1]
            double J15 = 1e-3 * (kf15 * IP3->Pcg * PIP2) / (Km15 / Cpc + PIP2); //[um^2-1 ms-1]
            double J16 = 1e-3 * kf16 * ec->IP3_state; //[uM * ms]

            //Update rate of change of state variables

            IP3R->Gd = (J7 + J13 + J12) - (J2 + J3);
            IP3R->Gt = J6 - (J7 + J9 + J10);
            IP3R->r = (-1.0) * (J1 + J2);
            IP3R->Rl = (J1 + J6) - J3;
            IP3R->Rg = J2 - J4;
            IP3R->Rlg = ((J3 - J5) + J4) - J6;
            IP3R->Rlgp = J5;
            IP3R->Pc = (J8 + J12) - J10;
            IP3R->Pcg = (J10 + J11) - J12;
            IP3R->P = J13 - (J9 + J8);
            IP3R->Pg = J9 - (J11 + J13);

            ecR->Ca_i += 1e-3 * Cpc * (-1.0) * (J8 + J11); // [uM/ms -> mM/ms] update Cai rate
            ecR->IP3_state = Cpc * (J14 + J15) - J16; //update IP3 rate [uM/ms]
        }
  }
}
