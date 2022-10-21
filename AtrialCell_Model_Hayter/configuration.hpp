#ifndef CONFIGURATION_H
#define CONFIGURATION_H

/*
Alter atrial_configuration to fit the required experiment.

All other cell parameters can be configured from within atrial_cell's member data, or within related ODE functions
*/

#include <string>

//Configuration class
struct atrial_configuration
{
    const std::string simulation_name = "RA"; //Name for simulation outputs
    const std::string initial_filename = "RA_1000_final_states.dat"; //File containing initial states to be parsed. Must be within "input_folder"
    const double dt = 0.001; //[ms] Time step for the simulation's solver
    const double output_dt = 0.01; //[ms] Intervals at which data is outputted. Must be divisible by dt
    const double temp = 310; //[K] Temperature for the simulation
    const bool take_measurements = true; //Toggle for AP measurement algorithm.

    //Configuration for steady stimulus
    const double amplitude = 10; //[pA/pF]
    const double duration = 4; //[ms]
    const double BCL = 1000; //[ms]
    const int beats = 10; //Number of heart beats for the simulation. Measurements require beats > 1

    //Channel blocks - Fraction of remaining current
    const double INa_Block = 1;
    const double ICaL_Block = 1;
    const double ICaT_Block = 1;
    const double IKur_Block = 1; //2.1 for LA
    const double Ito_Block = 1;
    const double IKr_Block = 1;
    const double IKACh_Block = 1;
    const double IKs_Block = 1;
    const double IKss_Block = 1;
    const double IK1_Block = 1; //1.7 for LA
    const double IKCa_Block = 1;
    const double IKb_Block = 1;
    const double If_Block = 1;
    const double INCX_Block = 1;
    const double ICab_Block = 1;
    const double INab_Block = 1;
    const double INaK_Block = 1;
    const double Jserca_Block = 1;
    const double JSRleak_Block = 1;

    //Configuration for other cell parameters:

    //IP3
    const double IP3_TIME = 0; //Time afterwhich IP3 has effect
    const bool IP3_dynamic_on = false;
    const double IP3_fixed = 0.023; //[uM] IP3 concentration - used if IP3 is not dynamic
    const double ET_1 = 0.1; //[uM] Ligand stimulation concentration - used if IP3 is dynamic. Test with 0.1 uM
    const bool IP3R_effect_on = false;

    //CaMKII
    const int _CaMKII_level_ = 0; //('WildType':0, 'OverExpression':1, '_Ko_ckOut':2)
    const int ACUTE = 0; //0:chronic OE, 1:acute OE (Does not require _CaMKII_level_ = 1 to have an effect)
    const double ACUTE_TIME = 1000; //Time afterwhich acute OE has effect

    //Caffeine
    const bool caffeine_on = false;
    const double caffeine_time = 16000; //[ms]

    //Isoproterenol
    const double _LigtotBA_ = 0; //[um] - SET LIGAND CONCENTRATION (0 or 0.1)
    const double _ISO_TIME_ = 1000; //&& t<BCL*100. Time afterwhich Isoproterenol has effect

    //ACh
    const double ACh = 0; //[uM]

    //Flags, with a value of 0 (false) or 1 (true)
    const int NaClampFlag = 0; //Na clamp => d_Nai_myo=0
    const int NaVsCaMKIIclamp = 0; //if 1, CaMKII Clamp on NaV
    const int PLMkoFlag = 0; //PLM KO (NaK 20% block, PLM_PKAp=100%)
    const int StrophFlag = 0; //Strophanthidin (dyastolic Na influx) => I_Nak=0
    const int DigitalisFlag = 0; //Digitalis (NaK 50% block)
    const int _LOOP_ = 0; //CaMKII-Na-Ca-CaMKII loop closed

    //Other configs
    const bool stim_with_K = true; //Stimulus current is potassium
    const bool fixed_Cl = false;
    const bool fixed_Ki = false;
    const bool serca_amp = false; //CaMKII effects on SERCA
    const bool faster_RyR = false;

    const bool new_INa = false;
    const bool new_INaK = false;

    const double max_m2_LCC_CK = 0.1;
};

#endif