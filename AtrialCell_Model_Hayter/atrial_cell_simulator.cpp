//Functions that operate the simulation

#include "atrial_cell.hpp"

//Prepares the simulation prior to solving the ODEs
void setup_simulation(atrial_cell &myocyte, const atrial_configuration &config)
{
    //Exception handling, ensuring valid user configuration
    try{
        //Setup simulation
        myocyte.output_config();
        myocyte.validate_config();
        myocyte.apply_config();
        myocyte.create_data_file();
        myocyte.create_measurement_file();
        myocyte.set_initial_states();
        std::cout << "Beginning " << config.simulation_name << " simulation..." << '\n' << std::endl;
    } catch(const std::exception &error) {
        std::cout << error.what() << std::endl; //Print error message associated with exception
        exit(EXIT_FAILURE);
    }
}

//Outputs simulation's configuration.
void atrial_cell::output_config()
{
    //Simulation info
    std::ofstream config_file{".\\output_folder\\" + config.simulation_name + "_config.txt"};
    if(config_file.fail()) throw std::runtime_error("output_folder not found in this directory"); //Check folder is present
    config_file << "Experiment type: " << "Steady Pacing" << '\n' << "\n";
    config_file << "Simulation name: " << config.simulation_name << '\n';
    config_file << "Initial states file name: " << config.initial_filename << '\n';
    config_file << "Final states file name: " << config.simulation_name << "_final_states.dat" << '\n' << "\n";

    config_file << "Solver time step(ms): " << config.dt << '\n';
    config_file << "Output time step(ms): " << config.output_dt << '\n';
    config_file << "Temperature(K): " << config.temp << '\n' << '\n';

    //Stimulus
    config_file << "Amplitude(mV): " << config.amplitude << '\n';
    config_file << "Duration(ms): " << config.duration << '\n';
    config_file << "BCL: " << config.BCL << '\n';
    config_file << "Frequency(Hz): " << freq << '\n';
    config_file << "Beats: " << config.beats << '\n' << "\n";

    //Channel blocks:
    const int number_of_blocks = 19;
    //Orders in the following two arrays must match up
    const std::string block_names[number_of_blocks] = {"INa Block", "ICaL Block", "ICaT Block", "IKur Block", "Ito Block", "IKr Block",
        "IKACh Block", "IKs Block", "IKss Block", "IK1 Block", "IKCa Block", "IKb Block", "INCX Block", "ICab Block", "INab Block",
        "INaK Block", "J SERCA Block", "J SR leak Block"};
    const double block_values[number_of_blocks] = {config.INa_Block, config.ICaL_Block, config.ICaT_Block, config.IKur_Block, config.Ito_Block,
        config.IKr_Block, config.IKACh_Block, config.IKs_Block, config.IKss_Block, config.IK1_Block, config.IKCa_Block, config.IKb_Block,
        config.INCX_Block, config.ICab_Block, config.INab_Block, config.INaK_Block, config.Jserca_Block, config.JSRleak_Block};
    //Ouput channel blocks that do not have the value of 1
    for(int i{}; i < number_of_blocks; i++) {
        if(i != 1) {
            config_file << block_names[i] << ": " << block_values[i] << '\n';
        }
    }

    //IP3:
    config_file << '\n' << "IP3_time: " << config.IP3_TIME << '\n';
    config_file << "Dynamic IP3: ";
    if(config.IP3_dynamic_on == false) { //IP3 fixed used
        config_file << "off" << '\n';
        config_file << "IP3 fixed(uM): " << config.IP3_fixed << '\n';
    } else { //ET_1 used
        config_file << "on" << '\n';
        config_file << "Ligand stimulation concentration (ET_1)(uM): " << config.ET_1 << '\n';
    }
    config_file << "IP3R effect: ";
    config.IP3R_effect_on == true ? config_file << "On" << '\n' : config_file << "off" << '\n' << '\n';

    //CaMKII
    config_file << "CaMKII level: ";
    if(config._CaMKII_level_ == 0) {
        config_file << "Wild type" << '\n';
    } else if(config._CaMKII_level_ == 1) {
        config_file << "Over expression" << '\n';
    } else if(config._CaMKII_level_ == 2) {
        config_file << "Knockout (KO)" << '\n';
    }
    config_file << "OE type: ";
    if(config.ACUTE == 0) {
        config_file << "Chronic" << '\n' << '\n';
    } else if(config.ACUTE == 1) {
        config_file << "Acute" << '\n' << '\n';
    }

    //Caffeine:
    config_file << "Caffeine effect: ";
    if(config.caffeine_on == true) {
        config_file << "on" << '\n';
        config_file << "Caffeine time(ms): " << config.caffeine_time << '\n' << '\n';
    } else {
        config_file << "off" << '\n' << '\n';
    }

    //Isoproterenol:
    if(config._LigtotBA_ == 0) {
        config_file << "Isoproterenol concentration: "<< config._LigtotBA_ << '\n' << '\n';
    } else {
        config_file << "Isoproterenol concentration: "<< config._LigtotBA_ << '\n';
        config_file << "ISO time: " << config._ISO_TIME_ << '\n' << '\n';
    }
    config_file << "ACh(uM): " << config.ACh << '\n' << '\n';

    //Flags
    config_file << "Flags:" << '\n';
    config_file << "Na Clamp: " << config.NaClampFlag << '\n';
    config_file << "CaMKII Clamp on NaV: " << config.NaVsCaMKIIclamp << '\n';
    config_file << "PLM KO: " << config.PLMkoFlag << '\n';
    config_file << "Strophanthidin: " << config.StrophFlag << '\n';
    config_file << "Digitalis: " << config.DigitalisFlag << '\n';
    config_file << "CaMKII-Na-Ca-CaMKII loop closed: " << config._LOOP_ << '\n' << '\n';

    //Other params:
    config_file << "Stim with K: ";
    config.stim_with_K == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Fixed Cl: ";
    config.fixed_Cl == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Fixed Ki: ";
    config.fixed_Ki == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "SERCA amp: ";
    config.fixed_Ki == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Faster RyR: ";
    config.serca_amp == true ? config_file << "on" << '\n' : config_file << "off" << '\n' << '\n';


    config_file << "New INa: ";
    config.new_INa == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "New INaK: ";
    config.new_INaK == true ? config_file << "on" << '\n' : config_file << "off" << '\n' << '\n';
    config_file << "Max m2 LCC CK: " << config.max_m2_LCC_CK << '\n' << '\n';
}

//Validates the user's core configuration setup
void atrial_cell::validate_config()
{
    //Check time-steps
    if(fabs(remainder(config.output_dt, config.dt)) > 1e-12) 
        throw std::logic_error("Invalid 'dt' and 'output_dt' configured. Check output_dt is divisible by dt.");
    if(config.output_dt <= 0 || config.dt <= 0) throw std::logic_error("Invalid 'dt' and 'output_dt' configured. Check that output_dt and dt are greater than 0.");
    //Check other params
    if(config.amplitude <= 0 || config.duration <=0 || config.BCL <=0 || config.beats <= 0) throw std::logic_error("Invalid steady pacing configuration. Ensure stimulus parameters are above zero");
}

//Applies some parameters dependent on user's configuration
void atrial_cell::apply_config()
{
    //Adjust CaMKII for mouse type
    if(config._CaMKII_level_ == 1) {
        _CaMKIItotJunc_ = 120.0*_n_OE_; //[um]
        _CaMKIItotSL_ = 120.0*8.293e-4*_n_OE_;
        _CaMKIItotCyt_ = 120.0*8.293e-4*_n_OE_;
    } else if(config._CaMKII_level_ == 2) {
        _CaMKIItotJunc_ = 0;
        _CaMKIItotSL_ = 0;
        _CaMKIItotCyt_ = 0;
    } else {
        _CaMKIItotJunc_ = 1 * 120.0; //[uM]
        _CaMKIItotSL_ = 1 * 120.0 * 8.293e-4;
        _CaMKIItotCyt_ = 1 * 120.0 * 8.293e-4;
    }

    if (config.new_INa) {
        _GNa_ *= 1.4;  // new I_Na (shanzhuo)
    } else {
        // _GNa_ *= 0.45; // not changed, to match experimental dvdtmax (weijian & shanzhuo)
    }

    if (config.faster_RyR) {
        _koCa_ *= 8;
        _kiCa_ *= 8;
        _kom_ *= 8;
        _kim_ *= 8;
    }
}

//Parse initial states datafile, validate and add to state vector
void atrial_cell::set_initial_states()
{
    double* state_ptr = &states.ec.V; //Pointer to first state

    std::ifstream initial_states_file;
    initial_states_file.open("./input_folder/" + config.initial_filename);
    if(initial_states_file.fail()) throw std::runtime_error("Initial states file not found, under ./input_folder/...");
    std::string line;
    //Loop for reading in states
    for(int i{}; i < states.number_of_states; i++) {
        if(getline(initial_states_file, line).fail() || (i + 1 == states.number_of_states && initial_states_file.peek() != -1)) { //Read in state and check number
            throw std::logic_error("Invalid initial states. " + std::to_string(i) + " states detected in initial states file but expected "
                + std::to_string(states.number_of_states) + ". Check initial states file or number_of_states variable definition.");
        }
        *state_ptr = atof(line.c_str()); //Convert string to double
        *(++state_ptr); //Move to next state
    }
    //Reset for IP3 module
    states.ip3.Gt = 0; //Gt
    states.ip3.Rl = 0; //Rl
    states.ip3.Rlg = 0; //Rlg
    states.ip3.Rlgp = 0; //Rlgp
    states.ip3.Pcg = 0; //Pcg
    states.ip3.Pg = 0; //Pg
}

//Determines whether stimulus is required and processes accordingly
void atrial_cell::process_stimulus(double t)
{
    double stim_temp;
    if(fmod(t, config.BCL) < config.duration){ //Condition for stimulus
        stim_temp = -config.amplitude;
    } else {
        stim_temp = 0; //Otherwise there is no stim for this t
    }
    if(abs(stim_temp - old_stim) > 1e-5){
        I_stim = stim_temp;
        old_stim = I_stim;
    }
}

//Updates the cells state vector via the forward Euler method
void atrial_cell::Euler_method()
{
    //For conciseness
    eccs *ec = &states.ec;
    bars *ba = &states.ba;
    ckiis *ck = &states.ck;
    cams *cj = &states.cj;
    cams *cs = &states.cs;
    cams *cc = &states.cc;
    ip3s *IP3 = &states.ip3;
    const eccs *ecR = &states_dot.ec;
    const bars *baR = &states_dot.ba;
    const ckiis *ckR = &states_dot.ck;
    const cams *cjR = &states_dot.cj;
    const cams *csR = &states_dot.cs;
    const cams *ccR = &states_dot.cc;
    const ip3s *IP3R = &states_dot.ip3;

    //Euler method
    ec->V += config.dt*ecR->V;
    ec->Nai_junc += config.dt*ecR->Nai_junc;
    ec->Nai_sl += config.dt*ecR->Nai_sl;
    ec->Nai_i += config.dt*ecR->Nai_i;
    ec->Ki += config.dt*ecR->Ki;
    ec->Ca_junc += config.dt*ecR->Ca_junc;
    ec->Ca_sl += config.dt*ecR->Ca_sl;
    ec->Ca_i += config.dt*ecR->Ca_i;
    ec->Ca_sr += config.dt*ecR->Ca_sr;
    ec->Cl_i += config.dt*ecR->Cl_i;
    ec->Csqn += config.dt*ecR->Csqn; 
    ec->m_na += config.dt*ecR->m_na; 
    ec->h_na += config.dt*ecR->h_na; 
    ec->j_na += config.dt*ecR->j_na; 
    ec->h_l += config.dt*ecR->h_l; 
    ec->R_SR += config.dt*ecR->R_SR; 
    ec->O_SR += config.dt*ecR->O_SR; 
    ec->I_SR += config.dt*ecR->I_SR; 
    ec->NaBj += config.dt*ecR->NaBj; 
    ec->NaBsl += config.dt*ecR->NaBsl; 
    ec->TnCL += config.dt*ecR->TnCL; 
    ec->TnCHc += config.dt*ecR->TnCHc; 
    ec->TnCHm += config.dt*ecR->TnCHm; 
    ec->Myosin_ca += config.dt*ecR->Myosin_ca; 
    ec->Myosin_mg += config.dt*ecR->Myosin_mg; 
    ec->SRB += config.dt*ecR->SRB; 
    ec->SLLj += config.dt*ecR->SLLj; 
    ec->SLLsl += config.dt*ecR->SLLsl; 
    ec->SLHj += config.dt*ecR->SLHj; 
    ec->SLHsl += config.dt*ecR->SLHsl; 
    ec->C2_m1j += config.dt*ecR->C2_m1j; 
    ec->C1_m1j += config.dt*ecR->C1_m1j; 
    ec->I1Ca_m1j += config.dt*ecR->I1Ca_m1j; 
    ec->I2Ca_m1j += config.dt*ecR->I2Ca_m1j; 
    ec->I1Ba_m1j += config.dt*ecR->I1Ba_m1j; 
    ec->I2Ba_m1j += config.dt*ecR->I2Ba_m1j; 
    ec->C2_m2j += config.dt*ecR->C2_m2j; 
    ec->C1_m2j += config.dt*ecR->C1_m2j; 
    ec->I1Ca_m2j += config.dt*ecR->I1Ca_m2j; 
    ec->I2Ca_m2j += config.dt*ecR->I2Ca_m2j; 
    ec->I1Ba_m2j += config.dt*ecR->I1Ba_m2j; 
    ec->I2Ba_m2j += config.dt*ecR->I2Ba_m2j; 
    ec->C2_m1sl += config.dt*ecR->C2_m1sl; 
    ec->C1_m1sl += config.dt*ecR->C1_m1sl; 
    ec->I1Ca_m1sl += config.dt*ecR->I1Ca_m1sl; 
    ec->I2Ca_m1sl += config.dt*ecR->I2Ca_m1sl; 
    ec->I1Ba_m1sl += config.dt*ecR->I1Ba_m1sl; 
    ec->I2Ba_m1sl += config.dt*ecR->I2Ba_m1sl; 
    ec->C2_m2sl += config.dt*ecR->C2_m2sl; 
    ec->C1_m2sl += config.dt*ecR->C1_m2sl; 
    ec->I1Ca_m2sl += config.dt*ecR->I1Ca_m2sl; 
    ec->I2Ca_m2sl += config.dt*ecR->I2Ca_m2sl; 
    ec->I1Ba_m2sl += config.dt*ecR->I1Ba_m2sl; 
    ec->I2Ba_m2sl += config.dt*ecR->I2Ba_m2sl; 
    ec->C_2ur += config.dt*ecR->C_2ur;
    ec->C_3ur += config.dt*ecR->C_3ur;
    ec->C_4ur += config.dt*ecR->C_4ur;
    ec->O_ur += config.dt*ecR->O_ur;
    ec->I_ur += config.dt*ecR->I_ur;
    ec->aKss += config.dt*ecR->aKss;
    ec->C3_Kr += config.dt*ecR->C3_Kr;
    ec->C1_Kr += config.dt*ecR->C1_Kr;
    ec->C2_Kr += config.dt*ecR->C2_Kr;
    ec->O_Kr += config.dt*ecR->O_Kr;
    ec->I_Kr_state += config.dt*ecR->I_Kr_state;
    ec->C2_to += config.dt*ecR->C2_to;
    ec->C1_to += config.dt*ecR->C1_to;
    ec->O_to += config.dt*ecR->O_to;
    ec->I_0to += config.dt*ecR->I_0to;
    ec->I_1to += config.dt*ecR->I_1to;
    ec->I_2to += config.dt*ecR->I_2to;
    ec->IP3r_Junc_O += config.dt*ecR->IP3r_Junc_O;
    ec->IP3r_Junc_RI += config.dt*ecR->IP3r_Junc_RI;
    ec->IP3r_Junc_R += config.dt*ecR->IP3r_Junc_R;
    ec->IP3r_Junc_RC += config.dt*ecR->IP3r_Junc_RC;
    ec->IP3r_Junc_RC += config.dt*ecR->IP3r_Junc_RC2;
    ec->IP3r_Junc_RC += config.dt*ecR->IP3r_Junc_RC3;
    ec->IP3_state += config.dt*ecR->IP3_state;
    ec->Ca_gap += config.dt*ecR->Ca_gap;
    ec->Ca_gapsr += config.dt*ecR->Ca_gapsr;
    ec->Csqn_gap += config.dt*ecR->Csqn_gap;
    ec->Ca_ss_gap += config.dt*ecR->Ca_ss_gap;
    ec->j_KACh += config.dt*ecR->j_KACh;
    ec->k_KACh += config.dt*ecR->k_KACh;
    ec->C2_KCa += config.dt*ecR->C2_KCa;
    ec->C3_KCa += config.dt*ecR->C3_KCa;
    ec->C4_KCa += config.dt*ecR->C4_KCa;
    ec->O1_KCa += config.dt*ecR->O1_KCa;
    ec->O2_KCa += config.dt*ecR->O2_KCa;
    cj->CaM += config.dt*cjR->CaM;
    cj->Ca2CaM += config.dt*cjR->Ca2CaM;
    cj->Ca4CaM += config.dt*cjR->Ca4CaM;
    cj->CaMB += config.dt*cjR->CaMB;
    cj->Ca2CaMB += config.dt*cjR->Ca2CaMB;
    cj->Ca4CaMB += config.dt*cjR->Ca4CaMB;
    cj->Pb2 += config.dt*cjR->Pb2;
    cj->Pb += config.dt*cjR->Pb;
    cj->Pt += config.dt*cjR->Pt;
    cj->Pt2 += config.dt*cjR->Pt2;
    cj->Pa += config.dt*cjR->Pa;
    cj->Ca4CaN += config.dt*cjR->Ca4CaN;
    cj->CaMCa4CaN += config.dt*cjR->CaMCa4CaN;
    cj->Ca2CaMCa4CaN += config.dt*cjR->Ca2CaMCa4CaN;
    cj->Ca4CaMCa4CaN += config.dt*cjR->Ca4CaMCa4CaN;
    cs->CaM += config.dt*csR->CaM;
    cs->Ca2CaM += config.dt*csR->Ca2CaM;
    cs->Ca4CaM += config.dt*csR->Ca4CaM;
    cs->CaMB += config.dt*csR->CaMB;
    cs->Ca2CaMB += config.dt*csR->Ca2CaMB;
    cs->Ca4CaMB += config.dt*csR->Ca4CaMB;
    cs->Pb2 += config.dt*csR->Pb2;
    cs->Pb += config.dt*csR->Pb;
    cs->Pt += config.dt*csR->Pt;
    cs->Pt2 += config.dt*csR->Pt2;
    cs->Pa += config.dt*csR->Pa;
    cs->Ca4CaN += config.dt*csR->Ca4CaN;
    cs->CaMCa4CaN += config.dt*csR->CaMCa4CaN;
    cs->Ca2CaMCa4CaN += config.dt*csR->Ca2CaMCa4CaN;
    cs->Ca4CaMCa4CaN += config.dt*csR->Ca4CaMCa4CaN;
    cc->CaM += config.dt*ccR->CaM;
    cc->Ca2CaM += config.dt*ccR->Ca2CaM;
    cc->Ca4CaM += config.dt*ccR->Ca4CaM;
    cc->CaMB += config.dt*ccR->CaMB;
    cc->Ca2CaMB += config.dt*ccR->Ca2CaMB;
    cc->Ca4CaMB += config.dt*ccR->Ca4CaMB;
    cc->Pb2 += config.dt*ccR->Pb2;
    cc->Pb += config.dt*ccR->Pb;
    cc->Pt += config.dt*ccR->Pt;
    cc->Pt2 += config.dt*ccR->Pt2;
    cc->Pa += config.dt*ccR->Pa;
    cc->Ca4CaN += config.dt*ccR->Ca4CaN;
    cc->CaMCa4CaN += config.dt*ccR->CaMCa4CaN;
    cc->Ca2CaMCa4CaN += config.dt*ccR->Ca2CaMCa4CaN;
    cc->Ca4CaMCa4CaN += config.dt*ccR->Ca4CaMCa4CaN;
    ck->LCC_CKjuncp += config.dt*ckR->LCC_CKjuncp ;
    ck->RyR2815p += config.dt*ckR->RyR2815p;
    ck->PLBT17p += config.dt*ckR->PLBT17p;
    ck->LCC_CKslp += config.dt*ckR->LCC_CKslp;
    ba->LR += config.dt*baR->LR;
    ba->LRG += config.dt*baR->LRG;
    ba->RG += config.dt*baR->RG;
    ba->b1AR_S464 += config.dt*baR->b1AR_S464;
    ba->b1AR_S301 += config.dt*baR->b1AR_S301;
    ba->GsaGTPtot += config.dt*baR->GsaGTPtot;
    ba->GsaGDP += config.dt*baR->GsaGDP;
    ba->Gsby += config.dt*baR->Gsby;
    ba->AC_GsaGTP += config.dt*baR->AC_GsaGTP;
    ba->PDEp += config.dt*baR->PDEp;
    ba->cAMPtot += config.dt*baR->cAMPtot;
    ba->RC_I += config.dt*baR->RC_I;
    ba->RCcAMP_I += config.dt*baR->RCcAMP_I;
    ba->RCcAMPcAMP_I += config.dt*baR->RCcAMPcAMP_I;
    ba->RcAMPcAMP_I += config.dt*baR->RcAMPcAMP_I;
    ba->PKACI += config.dt*baR->PKACI;
    ba->PKACI_PKI += config.dt*baR->PKACI_PKI;
    ba->RC_II += config.dt*baR->RC_II;
    ba->RCcAMP_II += config.dt*baR->RCcAMP_II;
    ba->RCcAMPcAMP_I += config.dt*baR->RCcAMPcAMP_II;
    ba->RcAMPcAMP_II += config.dt*baR->RcAMPcAMP_II;
    ba->PKACII += config.dt*baR->PKACII;
    ba->PKACII_PKI += config.dt*baR->PKACII_PKI;
    ba->I1p_PP1 += config.dt*baR->I1p_PP1;
    ba->I1ptot += config.dt*baR->I1ptot;
    ba->PLBp += config.dt*baR->PLBp;
    ba->PLMp += config.dt*baR->PLMp;
    ba->LCCap += config.dt*baR->LCCap;
    ba->LCCbp += config.dt*baR->LCCbp;
    ba->RyRp += config.dt*baR->RyRp;
    ba->TnIp += config.dt*baR->TnIp;
    ba->KURp += config.dt*baR->KURp;
    IP3->Gd += config.dt*IP3R->Gd;
    IP3->Gt += config.dt*IP3R->Gt;
    IP3->r += config.dt*IP3R->r;
    IP3->Rl += config.dt*IP3R->Rl;
    IP3->Rg += config.dt*IP3R->Rg;
    IP3->Rlg += config.dt*IP3R->Rlg;
    IP3->Rlgp += config.dt*IP3R->Rlgp;
    IP3->Pc += config.dt*IP3R->Pc;
    IP3->Pcg += config.dt*IP3R->Pcg;
    IP3->P += config.dt*IP3R->P;
    IP3->Pg += config.dt*IP3R->Pg;
}

//Take measurements at each simulation iteration.
//Measures APD start time from stimulus initiation.
void atrial_cell::measurements()
{
    //Progress to next beat after full BCL
    if(time == config.BCL*current_beat) { //Stabalised membrane potential, one step before stimulus
        RMP_prev = RMP; //Store beat's RMP before progressing to next beat (to be outputted)
        RMP = states.ec.V;
        APD_start = time; //APDs are measured from here
        output_meas = true; //Output beat's measurements
        take_measurements = true; //Start taking measurements for next beat
        current_beat +=1;
    }
    if(take_measurements == true && current_beat > 1 && config.take_measurements == true) { //Only measure useful data and skip first beat
        //Point at which max potential is reached within AP
        //Final condition ensures that the dV/dt transition is very likely to be within the AP
        if(dVdt_prev > 0 && dependents.dVdt <= 0 && time < config.BCL*(current_beat - 0.8)) { //Max potential at a change from positive to negative dV/dt
            max_V = V_prev;
            double upstroke_delta_V = max_V - RMP;
            APD90_V = upstroke_delta_V*0.1 + RMP;
            APD50_V = upstroke_delta_V*0.5 + RMP;
            APD30_V = upstroke_delta_V*0.7 + RMP;
        }
        if(states.ec.V <= APD30_V && V_prev > APD30_V) {
            APD30 = time - APD_start;
        }
        if(states.ec.V <= APD50_V && V_prev > APD50_V) {
            APD50 = time - APD_start;
        }
        if(states.ec.V <= APD90_V && V_prev > APD90_V) {
            APD90 = time - APD_start;
            norm_APD50 = APD50/APD90;
            take_measurements = false; //All APD measurements for given AP completed
        }
    }
}

//Calculate values required for each simulation loop
void atrial_cell::store_variables()
{
    dVdt_prev = dependents.dVdt;
    V_prev = states.ec.V;
}

void atrial_cell::update_time(const int &step)
{
    time = step * config.dt; //Update time
}

//Generates output file for core simulation data
void atrial_cell::create_data_file()
{
    //File to output currents, in a formatted string
    std::string file_path{".\\output_folder\\" + config.simulation_name + "_data.dat"};
    output_currents_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    //Output and format column labels
    fprintf(output_currents_file, "\t%-17s%-25s%-15s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-15s%-25s\n",
        "Time(ms)", "Membrane Potential(mV)", "dV/dt(V/s)", "Stimulus Current(pA/pF)", "I_Na(pA/pF)", "I_CaL(pA/pF)", "I_Ca_junc(pA/pF)",
        "I_Ca_sl(pA/pF)", "I_Kur(pA/pF)", "I_to(pA/pF)", "I_Kr(pA/pF)", "I_KACh(pA/pF)", "I_Kss(pA/pF)", "I_K1(pA/pF)", "I_KCa(pA/pF)",
        "I_Kb(pA/pF)", "I_NCX(pA/pF)", "I_Cap(pA/pF)", "I_NaK(pA/pF)", "I_Nab(pA/pF)", "I_Cab(pA/pF)",
        "I_ClCa(pA/pF)", "I_Clb(pA/pF)", "I_Na_fast(pA/pF)", "I_Na_slow(pA/pF)", "I_CaNa(pA/pF)", "I_CaK(pA/pF)", "J_SERCA(uM/ms)",
        "J_SRCarel(uM/ms)", "J_SRleak(uM/ms)", "Ca_junc(uM)", "Ca_sl(uM)", "Ca_i(uM)", "Ca_sr(uM)", "Na_Junc(mM)", "Na_sl(mM)", 
        "Nai_i(mM)", "K_i(mM)", "Cl_i(mM)", "CaMKII_junc_act(uM)", "CaMKII_sl_act(uM)", "CaMKII_cyt_act(uM)", "LCC_CKp(%)", "RyR_CKp(%)",
        "PLB_CKp_(%)", "PLM_PKAp(%)", "IP3(uM)", "IP3r_Junc_O");
}

//Generates output file for measurements on data
void atrial_cell::create_measurement_file()
{
    std::string file_path{".\\output_folder\\" + config.simulation_name + "_measurements.dat"};
    output_measurements_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    fprintf(output_measurements_file, "\t%-15s%-20s%-20s%-15s%-15s%-15s%-20s%-15s\n", "Beat Number", "AP Start Time(ms)", "Max Potential(mV)", "APD_30(ms)", "APD_50(ms)", "APD_90(ms)", "APD_50/APD_90", "Pre-AP RMP(mV)");
}

//Outputs data to relevant files on each simulation iteration. Use updated states, within the state vector
void atrial_cell::output_data()
{
    const eccs *ec = &states.ec;

    //Ensuring data is only outputted on the configured intervals
    if(fabs(remainder(time, config.output_dt)) > 1e-12) return;

    //If necessary, modify currents for output:

    const double I_CaNa_tot = dependents.I_CaNa_junc + dependents.I_CaNa_sl;
    const double I_CaK_tot = dependents.I_CaK_junc + dependents.I_CaK_sl;

    const double J_SERCA_scaled = dependents.J_serca*1e3; //[uM/s]
    const double J_SRCarel_scaled = dependents.J_SRCarel*1e3; //[uM/s]
    const double J_SRleak_scaled = dependents.J_SRleak*1e3; //[uM/s]

    const double Ca_junc_scaled = ec->Ca_junc*1e3; //[uM]
    const double Ca_sl_scaled = ec->Ca_sl*1e3; //[uM]
    const double Ca_i_scaled = ec->Ca_i*1e3; //[uM]
    const double Ca_sr_scaled = ec->Ca_sr*1e3; //[uM]

    const double LCC_CKp_perc = dependents.ph.LCC_CKp*100; //Convert from fractions to percentages
    const double RyR_CKp_perc = dependents.ph.RyR_CKp*100;
    const double PLB_CKp_perc = dependents.ph.PLB_CKp*100;
    const double PLM_PKAp_perc = dependents.ph.PLM_PKAp*100;

    //Output and format data
    //"-"" Left aligns. Before "."" gives character spacings between outputs
    //After "."" with a "g" gives the number of significant figures and displays in either standard form or 10^x form
    fprintf(output_currents_file, "\t%-17.10g%-25.10g%-15g%-25.8g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-15.10g%-25.10g\n",
        time, states.ec.V, dependents.dVdt, I_stim, dependents.I_Na, dependents.I_CaL, dependents.I_Ca_junc, dependents.I_Ca_sl,
        dependents.I_Kur, dependents.I_to, dependents.I_Kr, dependents.I_KACh, dependents.I_Kss, dependents.I_K1, dependents.I_KCa,
        dependents.I_Kb, dependents.I_NCX, dependents.I_Cap, dependents.I_NaK, dependents.I_Nab, dependents.I_Cab, dependents.I_ClCa,
        dependents.I_Clb, dependents.I_Na_fast, dependents.I_Na_slow, I_CaNa_tot, I_CaK_tot, J_SERCA_scaled, J_SRCarel_scaled,
        J_SRleak_scaled, Ca_junc_scaled, Ca_sl_scaled, Ca_i_scaled, Ca_sr_scaled, ec->Nai_junc, ec->Nai_sl, ec->Nai_i,
        ec->Ki, ec->Cl_i, dependents.camkii_junc_act, dependents.camkii_sl_act, dependents.camkii_cyt_act, LCC_CKp_perc,
        RyR_CKp_perc, PLB_CKp_perc, PLM_PKAp_perc, ec->IP3_state, ec->IP3r_Junc_O);
}

//Output measurements made throughout previous simulation beat
void atrial_cell::output_measurements()
{
    if(output_meas == true) { //Output after full beat
        std::cout << "Beat " << current_beat - 1 << " completed" << '\n';
        if(config.take_measurements == true) { //For configuration toggle
            if(current_beat > 2) { //Don't output measurements for first beat (1 beat ahead here)
                fprintf(output_measurements_file, "\t%-15i%-20.10g%-20.10g%-15.6g%-15.6g%-20.5g%-15.6g%-15.6g\n", current_beat - 1, APD_start, max_V, APD30, APD50, APD90, norm_APD50, RMP_prev);
                std::cout << "Max V: " << max_V << '\n';
                std::cout << "APD30: " << APD30 << '\n';
                std::cout << "APD50: " << APD50 << '\n';
                std::cout << "APD90: " << APD90 << '\n';
                std::cout << "APD50/APD90: " << norm_APD50 << '\n';
                std::cout << "Pre-AP RMP(mV): " << RMP_prev << '\n';
            }
        }
        std::cout << std::endl;
    }
    output_meas = false;
}

void atrial_cell::end_simulation()
{
    try{
        output_final_states();
    } catch(const std::exception &error) {
        std::cout << error.what() << std::endl; //Print error message associated with exception
        exit(EXIT_FAILURE);
    }
    std::cout << '\n' << "Simulation Completed";
}

//Output the final state variables
void atrial_cell::output_final_states()
{
    const double *state_ptr = &states.ec.V; //Pointer to the first state

    std::ofstream final_states_file{".\\output_folder\\" + config.simulation_name + "_final_states.dat"};
    if(final_states_file.fail()) throw std::runtime_error{"output_folder was not found at the end of the simulation."};
    final_states_file.precision(20);
    for(int i{}; i < states.number_of_states; i++) { //Output each state
        final_states_file << *state_ptr << std::endl;
        *(++state_ptr);
    }
}