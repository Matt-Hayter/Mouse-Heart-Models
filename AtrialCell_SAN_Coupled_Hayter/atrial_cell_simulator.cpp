//General functions that operate the atrial cell
#include "simulation_configuration.hpp"
#include "atrial_cell.hpp"

//Calls on all atrial cell ODEs
void atrial_cell::ODEs(const double &time)
{
    IP3_ODEs(time);
    cam_ODEs(time);
    camkii_ODEs(time);
    bars_ODEs(time);
    ecc_ODEs(time);
}

//Applies some parameters dependent on user's configuration
void atrial_cell::apply_config()
{
    //Adjust CaMKII for mouse type
    if(atrial_config._CaMKII_level_ == 1) {
        _CaMKIItotJunc_ = 120.0*_n_OE_; //[um]
        _CaMKIItotSL_ = 120.0*8.293e-4*_n_OE_;
        _CaMKIItotCyt_ = 120.0*8.293e-4*_n_OE_;
    } else if(atrial_config._CaMKII_level_ == 2) {
        _CaMKIItotJunc_ = 0;
        _CaMKIItotSL_ = 0;
        _CaMKIItotCyt_ = 0;
    } else {
        _CaMKIItotJunc_ = 1 * 120.0; //[uM]
        _CaMKIItotSL_ = 1 * 120.0 * 8.293e-4;
        _CaMKIItotCyt_ = 1 * 120.0 * 8.293e-4;
    }

    if (atrial_config.new_INa) {
        _GNa_ *= 1.4;  // new I_Na (shanzhuo)
    } else {
        // _GNa_ *= 0.45; // not changed, to match experimental dvdtmax (weijian & shanzhuo)
    }

    if (atrial_config.faster_RyR) {
        _koCa_ *= 8;
        _kiCa_ *= 8;
        _kom_ *= 8;
        _kim_ *= 8;
    }
}

void atrial_cell::output_config(std::ofstream &config_file)
{
    //Atrial cell configuration:
    config_file << "-----------------------Atrial Cell Configuration-----------------------" << '\n' << '\n';
    config_file << "Initial states file name: " << atrial_config.initial_filename << '\n';
    config_file << "Final states file name: " << sim_config.simulation_name << "_atrial_final_states.dat" << '\n' << "\n";
    //Channel blocks:
    const int number_of_blocks = 19;
    //Orders in the following two arrays must match up
    const std::string block_names[number_of_blocks] = {"INa Block", "ICaL Block", "ICaT Block", "IKur Block", "Ito Block", "IKr Block",
        "IKACh Block", "IKs Block", "IKss Block", "IK1 Block", "IKCa Block", "IKb Block", "If Block", "INCX Block", "ICab Block", "INab Block",
        "INaK Block", "J SERCA Block", "J SR leak Block"};
    const double block_values[number_of_blocks] = {atrial_config.INa_Block, atrial_config.ICaL_Block, atrial_config.ICaT_Block, atrial_config.IKur_Block, atrial_config.Ito_Block,
        atrial_config.IKr_Block, atrial_config.IKACh_Block, atrial_config.IKs_Block, atrial_config.IKss_Block, atrial_config.IK1_Block, atrial_config.IKCa_Block, atrial_config.IKb_Block,
        atrial_config.If_Block, atrial_config.INCX_Block, atrial_config.ICab_Block, atrial_config.INab_Block, atrial_config.INaK_Block, atrial_config.Jserca_Block, atrial_config.JSRleak_Block};
    //Ouput channel blocks that do not have the value of 1
    for(int i{}; i < number_of_blocks; i++) {
        if(i != 1) {
            config_file << block_names[i] << ": " << block_values[i] << '\n';
        }
    }

    //IP3:
    config_file << '\n' << "IP3_time: " << atrial_config.IP3_TIME << '\n';
    config_file << "Dynamic IP3: ";
    if(atrial_config.IP3_dynamic_on == false) { //IP3 fixed used
        config_file << "off" << '\n';
        config_file << "IP3 fixed(uM): " << atrial_config.IP3_fixed << '\n';
    } else { //ET_1 used
        config_file << "on" << '\n';
        config_file << "Ligand stimulation concentration (ET_1)(uM): " << atrial_config.ET_1 << '\n';
    }
    config_file << "IP3R effect: ";
    atrial_config.IP3R_effect_on == true ? config_file << "On" << '\n' : config_file << "off" << '\n' << '\n';

    //CaMKII
    config_file << "CaMKII level: ";
    if(atrial_config._CaMKII_level_ == 0) {
        config_file << "Wild type" << '\n';
    } else if(atrial_config._CaMKII_level_ == 1) {
        config_file << "Over expression" << '\n';
    } else if(atrial_config._CaMKII_level_ == 2) {
        config_file << "Knockout (KO)" << '\n';
    }
    config_file << "OE type: ";
    if(atrial_config.ACUTE == 0) {
        config_file << "Chronic" << '\n' << '\n';
    } else if(atrial_config.ACUTE == 1) {
        config_file << "Acute" << '\n' << '\n';
    }

    //Caffeine:
    config_file << "Caffeine effect: ";
    if(atrial_config.caffeine_on == true) {
        config_file << "on" << '\n';
        config_file << "Caffeine time(ms): " << atrial_config.caffeine_time << '\n' << '\n';
    } else {
        config_file << "off" << '\n' << '\n';
    }

    //Isoproterenol:
    if(atrial_config._LigtotBA_ == 0) {
        config_file << "Isoproterenol concentration: "<< atrial_config._LigtotBA_ << '\n' << '\n';
    } else {
        config_file << "Isoproterenol concentration: "<< atrial_config._LigtotBA_ << '\n';
        config_file << "ISO time: " << atrial_config._ISO_TIME_ << '\n' << '\n';
    }

    //Flags
    config_file << "Flags:" << '\n';
    config_file << "Na Clamp: " << atrial_config.NaClampFlag << '\n';
    config_file << "CaMKII Clamp on NaV: " << atrial_config.NaVsCaMKIIclamp << '\n';
    config_file << "PLM KO: " << atrial_config.PLMkoFlag << '\n';
    config_file << "Strophanthidin: " << atrial_config.StrophFlag << '\n';
    config_file << "Digitalis: " << atrial_config.DigitalisFlag << '\n';
    config_file << "CaMKII-Na-Ca-CaMKII loop closed: " << atrial_config._LOOP_ << '\n' << '\n';

    //Other params:
    config_file << "Stim with K: ";
    atrial_config.K_junction_current == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Fixed Cl: ";
    atrial_config.fixed_Cl == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Fixed Ki: ";
    atrial_config.fixed_Ki == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "SERCA amp: ";
    atrial_config.fixed_Ki == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "Faster RyR: ";
    atrial_config.serca_amp == true ? config_file << "on" << '\n' : config_file << "off" << '\n' << '\n';


    config_file << "New INa: ";
    atrial_config.new_INa == true ? config_file << "on" << '\n' : config_file << "off" << '\n';
    config_file << "New INaK: ";
    atrial_config.new_INaK == true ? config_file << "on" << '\n' : config_file << "off" << '\n' << '\n';
    config_file << "ACh(uM): " << atrial_config.ACh << '\n';
    config_file << "Max m2 LCC CK: " << atrial_config.max_m2_LCC_CK << '\n' << '\n';
}

//Parse initial states datafile, validate and add to state vector
void atrial_cell::set_initial_states()
{
    double* state_ptr = &states.ec.V; //Pointer to first state

    std::ifstream initial_states_file;
    initial_states_file.open("./input_folder/atrial_inputs/" + atrial_config.initial_filename);
    if(initial_states_file.fail()) throw std::runtime_error("Atrial initial states file not found, under ./input_folder/atrial_inputs/...");
    std::string line;
    //Loop for reading in states
    for(int i{}; i < states.number_of_states; i++) {
        if(getline(initial_states_file, line).fail() || (i + 1 == states.number_of_states && initial_states_file.peek() != -1)) { //Read in state and check number
            throw std::logic_error("Invalid atrial initial states. " + std::to_string(i) + " states detected in initial states file but expected "
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
    apply_config(); //Apply various configuration params
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

    //Backward Euler method
    ec->V += sim_config.dt*ecR->V;
    ec->Nai_junc += sim_config.dt*ecR->Nai_junc;
    ec->Nai_sl += sim_config.dt*ecR->Nai_sl;
    ec->Nai_i += sim_config.dt*ecR->Nai_i;
    ec->Ki += sim_config.dt*ecR->Ki;
    ec->Ca_junc += sim_config.dt*ecR->Ca_junc;
    ec->Ca_sl += sim_config.dt*ecR->Ca_sl;
    ec->Ca_i += sim_config.dt*ecR->Ca_i;
    ec->Ca_sr += sim_config.dt*ecR->Ca_sr;
    ec->Cl_i += sim_config.dt*ecR->Cl_i;
    ec->Csqn += sim_config.dt*ecR->Csqn; 
    ec->m_na += sim_config.dt*ecR->m_na; 
    ec->h_na += sim_config.dt*ecR->h_na; 
    ec->j_na += sim_config.dt*ecR->j_na; 
    ec->h_l += sim_config.dt*ecR->h_l; 
    ec->R_SR += sim_config.dt*ecR->R_SR; 
    ec->O_SR += sim_config.dt*ecR->O_SR; 
    ec->I_SR += sim_config.dt*ecR->I_SR; 
    ec->NaBj += sim_config.dt*ecR->NaBj; 
    ec->NaBsl += sim_config.dt*ecR->NaBsl; 
    ec->TnCL += sim_config.dt*ecR->TnCL; 
    ec->TnCHc += sim_config.dt*ecR->TnCHc; 
    ec->TnCHm += sim_config.dt*ecR->TnCHm; 
    ec->Myosin_ca += sim_config.dt*ecR->Myosin_ca; 
    ec->Myosin_mg += sim_config.dt*ecR->Myosin_mg; 
    ec->SRB += sim_config.dt*ecR->SRB; 
    ec->SLLj += sim_config.dt*ecR->SLLj; 
    ec->SLLsl += sim_config.dt*ecR->SLLsl; 
    ec->SLHj += sim_config.dt*ecR->SLHj; 
    ec->SLHsl += sim_config.dt*ecR->SLHsl; 
    ec->C2_m1j += sim_config.dt*ecR->C2_m1j; 
    ec->C1_m1j += sim_config.dt*ecR->C1_m1j; 
    ec->I1Ca_m1j += sim_config.dt*ecR->I1Ca_m1j; 
    ec->I2Ca_m1j += sim_config.dt*ecR->I2Ca_m1j; 
    ec->I1Ba_m1j += sim_config.dt*ecR->I1Ba_m1j; 
    ec->I2Ba_m1j += sim_config.dt*ecR->I2Ba_m1j; 
    ec->C2_m2j += sim_config.dt*ecR->C2_m2j; 
    ec->C1_m2j += sim_config.dt*ecR->C1_m2j; 
    ec->I1Ca_m2j += sim_config.dt*ecR->I1Ca_m2j; 
    ec->I2Ca_m2j += sim_config.dt*ecR->I2Ca_m2j; 
    ec->I1Ba_m2j += sim_config.dt*ecR->I1Ba_m2j; 
    ec->I2Ba_m2j += sim_config.dt*ecR->I2Ba_m2j; 
    ec->C2_m1sl += sim_config.dt*ecR->C2_m1sl; 
    ec->C1_m1sl += sim_config.dt*ecR->C1_m1sl; 
    ec->I1Ca_m1sl += sim_config.dt*ecR->I1Ca_m1sl; 
    ec->I2Ca_m1sl += sim_config.dt*ecR->I2Ca_m1sl; 
    ec->I1Ba_m1sl += sim_config.dt*ecR->I1Ba_m1sl; 
    ec->I2Ba_m1sl += sim_config.dt*ecR->I2Ba_m1sl; 
    ec->C2_m2sl += sim_config.dt*ecR->C2_m2sl; 
    ec->C1_m2sl += sim_config.dt*ecR->C1_m2sl; 
    ec->I1Ca_m2sl += sim_config.dt*ecR->I1Ca_m2sl; 
    ec->I2Ca_m2sl += sim_config.dt*ecR->I2Ca_m2sl; 
    ec->I1Ba_m2sl += sim_config.dt*ecR->I1Ba_m2sl; 
    ec->I2Ba_m2sl += sim_config.dt*ecR->I2Ba_m2sl; 
    ec->C_2ur += sim_config.dt*ecR->C_2ur;
    ec->C_3ur += sim_config.dt*ecR->C_3ur;
    ec->C_4ur += sim_config.dt*ecR->C_4ur;
    ec->O_ur += sim_config.dt*ecR->O_ur;
    ec->I_ur += sim_config.dt*ecR->I_ur;
    ec->aKss += sim_config.dt*ecR->aKss;
    ec->C3_Kr += sim_config.dt*ecR->C3_Kr;
    ec->C1_Kr += sim_config.dt*ecR->C1_Kr;
    ec->C2_Kr += sim_config.dt*ecR->C2_Kr;
    ec->O_Kr += sim_config.dt*ecR->O_Kr;
    ec->I_Kr_state += sim_config.dt*ecR->I_Kr_state;
    ec->C2_to += sim_config.dt*ecR->C2_to;
    ec->C1_to += sim_config.dt*ecR->C1_to;
    ec->O_to += sim_config.dt*ecR->O_to;
    ec->I_0to += sim_config.dt*ecR->I_0to;
    ec->I_1to += sim_config.dt*ecR->I_1to;
    ec->I_2to += sim_config.dt*ecR->I_2to;
    ec->IP3r_Junc_O += sim_config.dt*ecR->IP3r_Junc_O;
    ec->IP3r_Junc_RI += sim_config.dt*ecR->IP3r_Junc_RI;
    ec->IP3r_Junc_R += sim_config.dt*ecR->IP3r_Junc_R;
    ec->IP3r_Junc_RC += sim_config.dt*ecR->IP3r_Junc_RC;
    ec->IP3r_Junc_RC += sim_config.dt*ecR->IP3r_Junc_RC2;
    ec->IP3r_Junc_RC += sim_config.dt*ecR->IP3r_Junc_RC3;
    ec->IP3_state += sim_config.dt*ecR->IP3_state;
    ec->Ca_gap += sim_config.dt*ecR->Ca_gap;
    ec->Ca_gapsr += sim_config.dt*ecR->Ca_gapsr;
    ec->Csqn_gap += sim_config.dt*ecR->Csqn_gap;
    ec->Ca_ss_gap += sim_config.dt*ecR->Ca_ss_gap;
    ec->j_KACh += sim_config.dt*ecR->j_KACh;
    ec->k_KACh += sim_config.dt*ecR->k_KACh;
    ec->C2_KCa += sim_config.dt*ecR->C2_KCa;
    ec->C3_KCa += sim_config.dt*ecR->C3_KCa;
    ec->C4_KCa += sim_config.dt*ecR->C4_KCa;
    ec->O1_KCa += sim_config.dt*ecR->O1_KCa;
    ec->O2_KCa += sim_config.dt*ecR->O2_KCa;
    cj->CaM += sim_config.dt*cjR->CaM;
    cj->Ca2CaM += sim_config.dt*cjR->Ca2CaM;
    cj->Ca4CaM += sim_config.dt*cjR->Ca4CaM;
    cj->CaMB += sim_config.dt*cjR->CaMB;
    cj->Ca2CaMB += sim_config.dt*cjR->Ca2CaMB;
    cj->Ca4CaMB += sim_config.dt*cjR->Ca4CaMB;
    cj->Pb2 += sim_config.dt*cjR->Pb2;
    cj->Pb += sim_config.dt*cjR->Pb;
    cj->Pt += sim_config.dt*cjR->Pt;
    cj->Pt2 += sim_config.dt*cjR->Pt2;
    cj->Pa += sim_config.dt*cjR->Pa;
    cj->Ca4CaN += sim_config.dt*cjR->Ca4CaN;
    cj->CaMCa4CaN += sim_config.dt*cjR->CaMCa4CaN;
    cj->Ca2CaMCa4CaN += sim_config.dt*cjR->Ca2CaMCa4CaN;
    cj->Ca4CaMCa4CaN += sim_config.dt*cjR->Ca4CaMCa4CaN;
    cs->CaM += sim_config.dt*csR->CaM;
    cs->Ca2CaM += sim_config.dt*csR->Ca2CaM;
    cs->Ca4CaM += sim_config.dt*csR->Ca4CaM;
    cs->CaMB += sim_config.dt*csR->CaMB;
    cs->Ca2CaMB += sim_config.dt*csR->Ca2CaMB;
    cs->Ca4CaMB += sim_config.dt*csR->Ca4CaMB;
    cs->Pb2 += sim_config.dt*csR->Pb2;
    cs->Pb += sim_config.dt*csR->Pb;
    cs->Pt += sim_config.dt*csR->Pt;
    cs->Pt2 += sim_config.dt*csR->Pt2;
    cs->Pa += sim_config.dt*csR->Pa;
    cs->Ca4CaN += sim_config.dt*csR->Ca4CaN;
    cs->CaMCa4CaN += sim_config.dt*csR->CaMCa4CaN;
    cs->Ca2CaMCa4CaN += sim_config.dt*csR->Ca2CaMCa4CaN;
    cs->Ca4CaMCa4CaN += sim_config.dt*csR->Ca4CaMCa4CaN;
    cc->CaM += sim_config.dt*ccR->CaM;
    cc->Ca2CaM += sim_config.dt*ccR->Ca2CaM;
    cc->Ca4CaM += sim_config.dt*ccR->Ca4CaM;
    cc->CaMB += sim_config.dt*ccR->CaMB;
    cc->Ca2CaMB += sim_config.dt*ccR->Ca2CaMB;
    cc->Ca4CaMB += sim_config.dt*ccR->Ca4CaMB;
    cc->Pb2 += sim_config.dt*ccR->Pb2;
    cc->Pb += sim_config.dt*ccR->Pb;
    cc->Pt += sim_config.dt*ccR->Pt;
    cc->Pt2 += sim_config.dt*ccR->Pt2;
    cc->Pa += sim_config.dt*ccR->Pa;
    cc->Ca4CaN += sim_config.dt*ccR->Ca4CaN;
    cc->CaMCa4CaN += sim_config.dt*ccR->CaMCa4CaN;
    cc->Ca2CaMCa4CaN += sim_config.dt*ccR->Ca2CaMCa4CaN;
    cc->Ca4CaMCa4CaN += sim_config.dt*ccR->Ca4CaMCa4CaN;
    ck->LCC_CKjuncp += sim_config.dt*ckR->LCC_CKjuncp ;
    ck->RyR2815p += sim_config.dt*ckR->RyR2815p;
    ck->PLBT17p += sim_config.dt*ckR->PLBT17p;
    ck->LCC_CKslp += sim_config.dt*ckR->LCC_CKslp;
    ba->LR += sim_config.dt*baR->LR;
    ba->LRG += sim_config.dt*baR->LRG;
    ba->RG += sim_config.dt*baR->RG;
    ba->b1AR_S464 += sim_config.dt*baR->b1AR_S464;
    ba->b1AR_S301 += sim_config.dt*baR->b1AR_S301;
    ba->GsaGTPtot += sim_config.dt*baR->GsaGTPtot;
    ba->GsaGDP += sim_config.dt*baR->GsaGDP;
    ba->Gsby += sim_config.dt*baR->Gsby;
    ba->AC_GsaGTP += sim_config.dt*baR->AC_GsaGTP;
    ba->PDEp += sim_config.dt*baR->PDEp;
    ba->cAMPtot += sim_config.dt*baR->cAMPtot;
    ba->RC_I += sim_config.dt*baR->RC_I;
    ba->RCcAMP_I += sim_config.dt*baR->RCcAMP_I;
    ba->RCcAMPcAMP_I += sim_config.dt*baR->RCcAMPcAMP_I;
    ba->RcAMPcAMP_I += sim_config.dt*baR->RcAMPcAMP_I;
    ba->PKACI += sim_config.dt*baR->PKACI;
    ba->PKACI_PKI += sim_config.dt*baR->PKACI_PKI;
    ba->RC_II += sim_config.dt*baR->RC_II;
    ba->RCcAMP_II += sim_config.dt*baR->RCcAMP_II;
    ba->RCcAMPcAMP_I += sim_config.dt*baR->RCcAMPcAMP_II;
    ba->RcAMPcAMP_II += sim_config.dt*baR->RcAMPcAMP_II;
    ba->PKACII += sim_config.dt*baR->PKACII;
    ba->PKACII_PKI += sim_config.dt*baR->PKACII_PKI;
    ba->I1p_PP1 += sim_config.dt*baR->I1p_PP1;
    ba->I1ptot += sim_config.dt*baR->I1ptot;
    ba->PLBp += sim_config.dt*baR->PLBp;
    ba->PLMp += sim_config.dt*baR->PLMp;
    ba->LCCap += sim_config.dt*baR->LCCap;
    ba->LCCbp += sim_config.dt*baR->LCCbp;
    ba->RyRp += sim_config.dt*baR->RyRp;
    ba->TnIp += sim_config.dt*baR->TnIp;
    ba->KURp += sim_config.dt*baR->KURp;
    IP3->Gd += sim_config.dt*IP3R->Gd;
    IP3->Gt += sim_config.dt*IP3R->Gt;
    IP3->r += sim_config.dt*IP3R->r;
    IP3->Rl += sim_config.dt*IP3R->Rl;
    IP3->Rg += sim_config.dt*IP3R->Rg;
    IP3->Rlg += sim_config.dt*IP3R->Rlg;
    IP3->Rlgp += sim_config.dt*IP3R->Rlgp;
    IP3->Pc += sim_config.dt*IP3R->Pc;
    IP3->Pcg += sim_config.dt*IP3R->Pcg;
    IP3->P += sim_config.dt*IP3R->P;
    IP3->Pg += sim_config.dt*IP3R->Pg;
}

//Calculate values required for each simulation loop
void atrial_cell::store_variables()
{
    dVdt_old = dependents.dVdt; //Store data prior to Euler method
    V_old = states.ec.V;
}

//Take measurements at each simulation iteration.
//APDs now measured from point of max dV/dt
void atrial_cell::measurements(const double &time, int &current_beat)
{
    if(sim_config.take_measurements == true) { //If true, take measurements
        if(dVdt_old < 0 && dependents.dVdt >= 0) { //Measure MDP
            MDP = V_old;
            cycle_start = time;
            dVdt_max = -5000; //Reset for next iteration
        }
        if(dependents.dVdt>dVdt_max) {
            dVdt_max = dependents.dVdt;
            APD_start = time; //Start APD measurements at max dV/dt
        }
        if(current_beat > 1) { //First step of SAN AP and skip first beat
            //Point at which max potential is reached within AP
            //Final condition ensures that the dV/dt transition is very likely to be within the AP
            if(dVdt_old > 0 && dependents.dVdt <= 0) { //Max potential at a change from positive to negative dV/dt
                max_V = V_old;
                double upstroke_delta_V = max_V - MDP;
                APD90_V = upstroke_delta_V*0.1 + MDP;
                APD50_V = upstroke_delta_V*0.5 + MDP;
                APD30_V = upstroke_delta_V*0.7 + MDP;
            }
            if(states.ec.V <= APD30_V && V_old > APD30_V) {
                APD30 = time - APD_start;
            }
            if(states.ec.V <= APD50_V && V_old > APD50_V) {
                APD50 = time - APD_start;
            }
            if(states.ec.V <= APD90_V && V_old > APD90_V) {
                APD90 = time - APD_start;
                norm_APD50 = APD50/APD90;
                output_measurements(current_beat);
            }
        }
    }
}

//Output measurements made throughout previous simulation beat
void atrial_cell::output_measurements(const int &current_beat)
{
    fprintf(output_measurements_file, "\t%-15i%-20.10g%-20.6g%-15.6g%-15.5g%-20.6g%-15.6g%-15.6g\n", current_beat, cycle_start, max_V, MDP, dVdt_max, APD30, APD50, APD90, norm_APD50);
    std::cout << "\tAtrial Cell Measurements: " << '\n';
    std::cout << "\t\tMax potential(mV): " << max_V << std::endl;
    std::cout << "\t\tMDP(mV): " << MDP << std::endl;
    std::cout << "\t\tdV/dt max(V/s): " << dVdt_max << std::endl;
    std::cout << "\t\tAtrial APD30(ms): " << APD30 << '\n';
    std::cout << "\t\tAtrial APD50(ms): " << APD50 << '\n';
    std::cout << "\t\tAtrial APD90(ms): " << APD90 << '\n';
    std::cout << "\t\tAtrial APD50/APD90: " << norm_APD50 << '\n';
}

//Generates output file for core simulation data
void atrial_cell::create_data_file()
{
    //File to output currents, in a formatted string
    std::string file_path{".\\output_folder\\" + sim_config.simulation_name + "_atrial_data.dat"};
    output_currents_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    //Output and format column labels
    fprintf(output_currents_file, "\t%-17s%-25s%-15s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-15s%-25s%-25s\n",
        "Time(ms)", "Membrane Potential(mV)", "dV/dt(V/s)", "Junction Current(pA/pF)", "I_Na(pA/pF)", "I_CaL(pA/pF)", "I_Ca_junc(pA/pF)",
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
    std::string file_path{".\\output_folder\\" + sim_config.simulation_name + "_atrial_measurements.dat"};
    output_measurements_file = fopen(file_path.c_str(), "w");
    if(output_currents_file == nullptr) throw std::runtime_error("output_folder not found in this directory"); //Check output folder is present
    fprintf(output_measurements_file, "\t%-15s%-20s%-20s%-15s%-15s%-20s%-15s%-15s\n", "Beat Number", "AP Start Time(ms)", "Max Potential(mV)", "MDP(mV)", "Max dV/dt(V/s)", "APD_30(ms)", "APD_50(ms)", "APD_90(ms)", "APD_50/APD_90");
}

//Outputs data to relevant files on each simulation iteration. Use updated states, within the state vector
void atrial_cell::output_data(const double &time)
{
    const eccs *ec = &states.ec;

    //Ensuring data is only outputted on the configured intervals
    if(fabs(remainder(time, sim_config.output_dt)) > 1e-12) return;

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
    fprintf(output_currents_file, "\t%-17.10g%-25.10g%-15g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-25.10g%-15.10g%-25.10g%-25.10g\n",
        time, states.ec.V, dependents.dVdt, I_j, dependents.I_Na, dependents.I_CaL, dependents.I_Ca_junc, dependents.I_Ca_sl,
        dependents.I_Kur, dependents.I_to, dependents.I_Kr, dependents.I_KACh, dependents.I_Kss, dependents.I_K1, dependents.I_KCa,
        dependents.I_Kb, dependents.I_NCX, dependents.I_Cap, dependents.I_NaK, dependents.I_Nab, dependents.I_Cab, dependents.I_ClCa,
        dependents.I_Clb, dependents.I_Na_fast, dependents.I_Na_slow, I_CaNa_tot, I_CaK_tot, J_SERCA_scaled, J_SRCarel_scaled,
        J_SRleak_scaled, Ca_junc_scaled, Ca_sl_scaled, Ca_i_scaled, Ca_sr_scaled, ec->Nai_junc, ec->Nai_sl, ec->Nai_i,
        ec->Ki, ec->Cl_i, dependents.camkii_junc_act, dependents.camkii_sl_act, dependents.camkii_cyt_act, LCC_CKp_perc,
        RyR_CKp_perc, PLB_CKp_perc, PLM_PKAp_perc, ec->IP3_state, ec->IP3r_Junc_O);
}

//Output the final state variables
void atrial_cell::output_final_states()
{
    const double *state_ptr = &states.ec.V; //Pointer to the first state

    std::ofstream final_states_file{".\\output_folder\\" + sim_config.simulation_name + "_atrial_final_states.dat"};
    if(final_states_file.fail()) throw std::runtime_error{"output_folder was not found at the end of the simulation."};
    final_states_file.precision(20);
    for(int i{}; i < states.number_of_states; i++) { //Output each state
        final_states_file << *state_ptr << std::endl;
        *(++state_ptr);
    }
}