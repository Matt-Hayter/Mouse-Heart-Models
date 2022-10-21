#include "atrial_cell.hpp"

void atrial_cell::ecc_ODEs()
{
    const eccs *ec = &states.ec;
    eccs *ecR = &states_dot.ec;

    dependents.GNaB = _GNaB_;
    dependents.Ena_junc = (1 / FoRT) * log(_Nao_ / ec->Nai_junc);            // [mV]
    dependents.Ena_sl = (1 / FoRT) * log(_Nao_ / ec->Nai_sl);            // [mV]
    dependents.Ek = (1 / FoRT) * log(_Ko_ / ec->Ki);            // [mV]
    dependents.Eca_junc = (1 / FoRT / 2) * log(_Cao_ / ec->Ca_junc);            // [mV]
    dependents.Eca_sl = (1 / FoRT / 2) * log(_Cao_ / ec->Ca_sl);            // [mV]
    dependents.Ecl = (1 / FoRT) * log(ec->Cl_i / _Clo_);            // [mV]

    //Calculate ecc currents
    double I_Na = INa_channel();
    double I_Nak = INaK_channel();
    double I_to = config.Ito_Block * Ito_channel();
    double I_Kur = config.IKur_Block * IKur_channel();
    double I_Kss = config.IKss_Block * IKss_channel();
    double I_Kr = config.IKr_Block * IKr_channel();
    double I_KACh = config. IKACh_Block * IKACh_channel();
    double I_K1 = config.IK1_Block * IK1_channel();
    double I_KCa = config.IKCa_Block * IKCa_channel();
    double I_Kb = config.IKb_Block * background_currents(6);
    double I_NCX = INCX_channel();
    double I_Nab = background_currents(1);
    double I_ClCa = background_currents(2);
    double I_Clb = background_currents(3);
    double I_Cap = background_currents(4);
    double I_Cab = background_currents(5);
    double I_CaL = LCC_channel();

    Na_concentrations();
    classic_Ca_handling();

    //// Total Currents [uA/uF]
    double I_Na_tot_junc = dependents.I_Na_junc + dependents.I_Nab_junc + 3 * dependents.I_NCX_junc
        + 3 * dependents.I_Nak_junc + dependents.I_CaNa_junc;
    double I_Na_tot_sl = dependents.I_Na_sl + dependents.I_Nab_sl + 3 * dependents.I_NCX_sl
        + 3 * dependents.I_Nak_sl + dependents.I_CaNa_sl;

    double I_Ca_tot_junc =dependents.I_Ca_junc + dependents.I_Cab_junc + dependents.I_Cap_junc
        - 2 * dependents.I_NCX_junc;    // [uA/uF]
    double I_Ca_tot_sl = dependents.I_Ca_sl + dependents.I_Cab_sl + dependents.I_Cap_sl
        - 2 * dependents.I_NCX_sl;    // [uA/uF]

    double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
    double I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
    double I_K_tot = I_to + I_Kr + I_KACh + I_K1 + I_KCa + I_Kb - 2 * I_Nak + I_Kur + I_Kss + dependents.I_CaK_junc
        + dependents.I_CaK_sl;  //[mM/msec]

    double I_Cl_tot = I_ClCa + I_Clb;

    double I_tot = I_Na_tot + I_Ca_tot + I_K_tot + I_Cl_tot;

    // Update membrane voltage rate
    ecR->V = -(I_tot + I_stim);

    //Control of intracellular [K]
    if(config.fixed_Ki) {
        ecR->Ki = 0;
    } else {
        //intracellular [K] set to fixed or updated by this formula [mM/msec]
        if (config.stim_with_K) {
            ecR->Ki = -(I_K_tot + I_stim) * Cmem / (Vmyo * Frdy);
        } else {
            ecR->Ki = -(I_K_tot) * Cmem / (Vmyo * Frdy);
        }
    }

    //Control of intraceullular [Cl]
    if(config.fixed_Cl) {
        ecR->Cl_i = 0;
    } else {
        ecR->Cl_i = I_Cl_tot * Cmem / (Vmyo * Frdy);

    }

    //Update ecc currents to be outputted, along with
    dependents.dVdt = ecR->V;
    dependents.I_Na = I_Na;
    dependents.I_NaK = I_Nak;
    dependents.I_to = I_to;
    dependents.I_Kur = I_Kur;
    dependents.I_Kss = I_Kss;
    dependents.I_Kr = I_Kr;
    dependents.I_KACh = I_KACh;
    dependents.I_K1 = I_K1;
    dependents.I_KCa = I_KCa;
    dependents.I_Kb = I_Kb;
    dependents.I_NCX = I_NCX;
    dependents.I_Nab = I_Nab;
    dependents.I_ClCa = I_ClCa;
    dependents.I_Clb = I_Clb;
    dependents.I_Cap = I_Cap;
    dependents.I_Cab = I_Cab;
    dependents.I_CaL = I_CaL;

}

//Compute current from the INa channel
double atrial_cell::INa_channel()
{
    const auto &ec = states.ec;
    auto &ecR = states_dot.ec;
    double V = ec.V;
    double GNa = _GNa_; // [mS/uF]

    // Max INa alterations with CaMKII hyperactivity as in Hund & Rudy 2008
    double inashift = 0;
    double alphaCKII = 0;
    double deltGbarNal_CKII = 0; // no Na Gain in OE

    if (config._CaMKII_level_ == 1 && (config.ACUTE == 1 && time > config.ACUTE_TIME)) // acute CaMKII effects
    {
        inashift = -3.25;
        alphaCKII = -0.18;
        deltGbarNal_CKII = 3;
    }
    if (config._LOOP_ == 1) {
        double RyRp_WT_mean = 0.2101;
        double RyRp_OE_mean = 0.7387;      // Derived (1 Hz, no loop)
        double RyRp_OEloop_min = 0.7033;   // Derived (1 Hz, OE loop)
        double delta_loop = (3.0 / (RyRp_OE_mean - RyRp_WT_mean)) * dependents.ph.RyR_CKp - (3.0 / (RyRp_OE_mean - RyRp_WT_mean)) * RyRp_WT_mean;

        if (config.NaVsCaMKIIclamp == 1) {
        delta_loop = (3.0 / (RyRp_OE_mean - RyRp_WT_mean)) * RyRp_OEloop_min - (3.0 / (RyRp_OE_mean - RyRp_WT_mean)) * RyRp_WT_mean;
        }
        if (config._CaMKII_level_ == 1 && (config.ACUTE == 1 && time > config.ACUTE_TIME)) deltGbarNal_CKII = delta_loop; //acute CaMKII effects
        dependents.GNaB = _GNaB_ * (1 + delta_loop);
    }
    double I_Na_fast_junc;
    double I_Na_fast_sl;
    double m_inf, h_inf, j_inf;
    double tau_m, tau_h, tau_j;
    if (config.new_INa) {
    V = V - inashift;
    ///////////////////// New I_Na ////////////////////
    // m and h infinity changed according to data in shekhar 2018
    if (config.temp == 310) {
        // fit to temperature shifted data from 22C to 37C
        // m_inf, h_inf, tau_m totally refitted, the right branch of tau_h shifted with a Q10 factor
        // m_inf = 1 / (1 + exp(-(V + 49.8) / 8.2));
        // h_inf = 1 / (1 + exp((V + 74.5) / 8.3));
        // j_inf = 1 / (1 + exp((V + 60.1) / 3.3));
        m_inf = 1 / (1 + exp(-(V + 49) / 7.5));
        h_inf = 1 / (1 + exp((V + 76.6) / 6.8));
        // if (V > -55) h_inf = 1 / (1 + exp((V + 70.1) / 4.8));
        j_inf = h_inf;
        tau_m = (1 / (5 * exp((V + 25) / 27) + 13 / (1 + exp((-V - 45) / 7)) + 8.6 * exp(-(V + 60.4) / 6)));
        double left = 1.2 * exp(-(V + 130.2) / 15.3);
        double right = 0.7 * exp((V + 39.5) / 15.27) / pow(Q10_tau_h_Na, Qpow);
        // tau_h = 0.7 / (left + right) + 0.35;
        tau_h = 1.4 / (left + right) + 0.35;
        double left_j = 2 * 1.2 * exp(-(V + 72.2) / 3.3);
        double right_j = 1 * 0.4 * exp((V + 34.5) / 4.3);
        tau_j = 1.3 / (left_j + right_j) + 1;
        // double left_j = 0.8 * exp(-(V + 80.2) / 7.3);
        // double right_j = 1 * 0.4 * exp((V + 31.5) / 7.3);
        // tau_j = 1.05 / (left_j + right_j) + 2;
    } else {
        // fit to data at room temperature
        m_inf = 1 / (1 + exp(-(V + 60.57) / 8.5));
        h_inf = 1 / (1 + exp((V + 90) / 7.5));
        j_inf = h_inf;
        tau_m = 6 * 1 / (1 * exp((V + 75.6) / 10.8) + 8.6 * exp(-(V + 57.4) / 3)) + 0.1;
        tau_h = 2 * 1 / (1.2 * exp(-(V + 130.2) / 15.3) + exp((V + 39.5) / 15.3));
        tau_j = 2 * 1 / (1.2 * exp(-(V + 130.2) / 15.3) + exp((V + 39.5) / 15.3));
    }
    ecR.m_na = (m_inf - ec.m_na) / tau_m;
    ecR.h_na = (h_inf - ec.h_na) / tau_h;
    ecR.j_na = (j_inf - ec.j_na) / tau_j;
    I_Na_fast_junc = config.INa_Block * Fjunc_na * GNa * pow(ec.m_na, 3) * ec.h_na * ec.j_na * (V - dependents.Ena_junc);
    I_Na_fast_sl = config.INa_Block * Fsl_na * GNa * pow(ec.m_na, 3) * ec.h_na * ec.j_na * (V - dependents.Ena_sl);
    } else {
        /////////////////////// old I_Na /////////////////////////////
        double am = 0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
        double bm = 0.08 * exp(-(V) / 11);
        double ah = 0.135 * exp((80 + (V - inashift)) / -6.8);
        double bh = 1.1 * 3.56 * exp(0.079 * (V - inashift - 2)) + 3.1e5 * exp(0.35 * (V - inashift - 2));
        double aj = (1 + alphaCKII) * ((-1.2714e5 * exp(0.2444 * (V - inashift)) - 3.474e-5 * exp(-0.04391 * (V - inashift)))
            * ((V - inashift) + 37.78)  / (1 + exp(0.311 * ((V - inashift) + 79.23))));
        double bj = 0.1212 * exp(-0.01052 * (V - inashift)) / (1 + exp(-0.1378 * ((V - inashift) + 40.14)));
    if ((V - inashift) >= -40) {
        ah = 0;
        aj = 0;

        bh = 0.66 * 1 / (0.13 * (1 + exp(-((V - inashift) + 10.66) / 11.1))); // MOUSE
        bj = 0.3 * exp(-2.535e-7 * (V - inashift)) / (1 + exp(-0.1 * ((V - inashift) + 32)));
    }
    if ((V - inashift) < -40) {
        aj *= 5;
        bj *= 5;
    }
    ecR.m_na = (am * (1 - ec.m_na) - bm * ec.m_na);
    ecR.h_na = (ah * (1 - ec.h_na) - bh * ec.h_na);
    ecR.j_na = (aj * (1 - ec.j_na) - bj * ec.j_na);

    I_Na_fast_junc = config.INa_Block * Fjunc_na * GNa * pow(ec.m_na, 3) * ec.h_na * ec.j_na * (V - dependents.Ena_junc);
    I_Na_fast_sl = config.INa_Block * Fsl_na * GNa * pow(ec.m_na, 3) * ec.h_na * ec.j_na * (V - dependents.Ena_sl);
    }

    //I_Na,L: Late INa current (as in Hund & Rudy 2008)
    double GbarNal = 2 * 0.0065 * (1 + deltGbarNal_CKII); // deltGbar assigned in 'Fast INa' section

    double hlss = 1 / (1 + exp((V + 91) / 6.1));

    if (config._CaMKII_level_ == 1 ||
        (config.ACUTE == 1 && time > config.ACUTE_TIME)) //weijian, acute and chronic CaMKII effects
        hlss = 1 / (1 + exp((V + 91 - 6.8) / 6.1));

    double tauhl = 600.0; // ms
    ecR.h_l = (hlss - ec.h_l) / tauhl;
    double I_Na_slow_junc = config.INa_Block * Fjunc_na * GbarNal * pow(ec.m_na, 3) * ec.h_l * (V - dependents.Ena_junc);
    double I_Na_slow_sl = config.INa_Block * Fsl_na * GbarNal * pow(ec.m_na, 3) * ec.h_l * (V - dependents.Ena_sl);

    //I_Na: compute total current (fast and late components)
    double I_Na_junc = I_Na_fast_junc + I_Na_slow_junc;
    double I_Na_sl = I_Na_fast_sl + I_Na_slow_sl;

    dependents.I_Na_junc = I_Na_junc;
    dependents.I_Na_sl = I_Na_sl;
    dependents.I_Na_fast = I_Na_fast_junc + I_Na_fast_sl;
    dependents.I_Na_slow = I_Na_slow_junc + I_Na_slow_sl;

    return (I_Na_junc + I_Na_sl); //Return total Na currents
}

//Compute current from the INaK channel
double atrial_cell::INaK_channel()
{
    const eccs *ec = &(states.ec);
    double V = ec->V;

    double IbarNaK = _IbarNaK_; // [uA/uF] changed from rabbit (1.90719)
    double PLM_PKAp = dependents.ph.PLM_PKAp;

    if (config._CaMKII_level_ == 1 || (config.ACUTE == 1 && time > config.ACUTE_TIME)) IbarNaK = IbarNaK * 0.9; // acute and chronic CaMKII effects
    if (config.DigitalisFlag == 1) IbarNaK = IbarNaK * 0.5; // 50% block
    if (config.PLMkoFlag == 1) {
        PLM_PKAp = 1;
        IbarNaK = IbarNaK * 0.8;
    }
    if (config.StrophFlag == 1) IbarNaK = 0;

    double I_Nak;
    if (!config.new_INaK) {
        double _KmNai_p = 19.0; // [mM] !!! changed from 19 (morotti mouse ventricle) to 11
        double Km_Ko_ = 1.5; // [mM]

        double sigma = (exp(_Nao_ / 67.3) - 1) / 7.0;
        double fnak = 1 / (1 + 0.1245 * exp(-0.1 * V * FoRT) + 0.0365 * sigma * exp(-V * FoRT));
        double fracPKA_PLMo = 0.116738; // Derived quantity (PLM_PKAp(baseline)/PLMtot)
        double fracPKA_PLMiso = 0.859251; // Derived quantity (PLM_PKAp(ISO)/PLMtot)
        double kPKA_PLM = _KmNai_p * (1 - 0.7019) / (fracPKA_PLMiso / fracPKA_PLMo - 1); // PLM_PKAp ISO
        double _KmNai_p_PKA = -kPKA_PLM + kPKA_PLM * (PLM_PKAp / fracPKA_PLMo);
        _KmNai_p = _KmNai_p - _KmNai_p_PKA;

        dependents.I_Nak_junc = config.INaK_Block * Fjunc_nak * IbarNaK * fnak * _Ko_ / (1 + pow((_KmNai_p / ec->Nai_junc), 4)) / (_Ko_ + Km_Ko_);
        dependents.I_Nak_sl = config.INaK_Block * Fsl_nak * IbarNaK * fnak * _Ko_ / (1 + pow((_KmNai_p / ec->Nai_sl), 4)) / (_Ko_ + Km_Ko_);
        I_Nak = (dependents.I_Nak_junc + dependents.I_Nak_sl);

    } else {
        // I_NaK in the sub-sarcolemma area
        double v = ec->V;
        double nai = ec->Nai_sl;
        double ki = ec->Ki;

        double k1p = 949.5;
        double k1m = 182.4;
        double k2p = 687.2;

        double k2m = 39.4;
        double k3p = 1899.0;
        double k3m = 79300.0;
        double k4p = 639.0;
        double k4m = 40.0;
        double Knai0 = 9.073;
        double Knao0 = 27.78;
        double delta = -0.1550;
        double Knai = Knai0 * exp(FoRT * (delta * v) / 3.0);
        double Knao = Knao0 * exp(FoRT * ((1.0 - delta) * v) / 3.0);
        double Kki = 0.5;
        double Kko = 0.3582;
        double MgADP = 0.05;
        double MgATP = 9.8;
        double Kmgatp = 1.698e-7;
        double H = 1.0e-7;
        double eP = 4.2;
        double Khp = 1.698e-7;
        double Knap = 224.0;
        double Kxkur = 292.0;

        double P = eP / (1.0 + H / Khp + nai / Knap + ki / Kxkur);
        double a1 = (k1p * pow(nai / Knai, 3.0)) / (pow(1.0 + nai / Knai, 3.0) + pow(1.0 + ki / Kki, 2.0) - 1.0);
        double b1 = k1m * MgADP;
        double a2 = k2p;
        double b2 = (k2m * pow(_Nao_ / Knao, 3.0)) / (pow(1.0 + _Nao_ / Knao, 3.0) + pow(1.0 + _Ko_ / Kko, 2.0) - 1.0);
        double a3 = (k3p * pow(_Ko_ / Kko, 2.0)) / (pow(1.0 + _Nao_ / Knao, 3.0) + pow(1.0 + _Ko_ / Kko, 2.0) - 1.0);
        double b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
        double a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
        double b4 = (k4m * pow(ki / Kki, 2.0)) / (pow(1.0 + nai / Knai, 3.0) + pow(1.0 + ki / Kki, 2.0) - 1.0);

        double x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
        double x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
        double x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
        double x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;

        double E1 = x1 / (x1 + x2 + x3 + x4);
        double E2 = x2 / (x1 + x2 + x3 + x4);
        double E3 = x3 / (x1 + x2 + x3 + x4);
        double E4 = x4 / (x1 + x2 + x3 + x4);

        double zk = 1.0, zna = 1.0;
        double JnakNa = 3.0 * (E1 * a3 - E2 * b3);
        double JnakK = 2.0 * (E4 * b1 - E3 * a1);

        double Pnak = 30;
        I_Nak = config.INaK_Block * Pnak * (zna * JnakNa + zk * JnakK);

        dependents.I_Nak_junc = 0;
        dependents.I_Nak_sl = I_Nak;
    }
    return I_Nak;
}

//Compute current from the Ito channel
double atrial_cell::Ito_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;
    double I_to;

    double GtoFast = 0.44;
    
    //I_to MM Model
    double C_3to = 1 - ec->C2_to - ec->C1_to - ec->I_2to - ec->I_1to - ec->I_0to - ec->O_to;

    double _act = 1.0; //pow(5.4, 1.5);
    double _inact = 1.0; //pow(1.9, 1.5);
    double _iv = 1.0; //pow(1.37, 1.5);

    //Update transition rates  used to calculate ion channel current states in net step
    double a_to = 0.250270 * (exp(0.787791 * (V + 5) * FoRT)) * KQ10_to * _act;
    double b_to = 0.168346 * (exp(-0.757095 * (V + 46) * FoRT)) * KQ10_to * _act;

    double K_fto = 0.0233856 * KQ10_to * _inact;
    double K_bto = (K_fto / 40.0) * KQ10_to * _inact;

    double K_f2to = 0.430214 * KQ10_to * _inact;
    double K_b2to = ((360.000 * K_f2to * K_bto) / K_fto) * KQ10_to * _inact;

    I_to = pow(Q10to_conductance, Qpow) * g_Kto * ec->O_to * (V - dependents.Ek) * _iv;

    if (config._CaMKII_level_ == 1 || (config.ACUTE == 1 && time > config.ACUTE_TIME)) K_fto = K_fto * 5; // chronic and acute CaMKII effects

    // I_to=I_to*2.0/3.0; (chronic effects)

    //Update rates for probabilities of being in each Markov state, depending on transition rates
    ecR->C2_to = ((((3 * a_to * C_3to - b_to * ec->C2_to) + 2 * b_to * ec->C1_to) - 2 * a_to * ec->C2_to)
        + K_b2to * ec->I_2to) - K_f2to * ec->C2_to;
    ecR->C1_to = ((((2 * a_to * ec->C2_to - 2 * b_to * ec->C1_to) + 3 * b_to * ec->O_to) - a_to * ec->C1_to)
        + K_bto * ec->I_1to) - K_fto * ec->C1_to;
    ecR->O_to = ((a_to * ec->C1_to - 3 * b_to * ec->O_to) + K_bto * ec->I_0to) - K_fto * ec->O_to;
    ecR->I_0to = ((a_to * ec->I_1to - 3 * b_to * ec->I_0to) + K_fto * ec->O_to) - K_bto * ec->I_0to;
    ecR->I_1to = ((((2 * a_to * ec->I_2to - (2 * b_to * ec->I_1to) / 360) + 3 * b_to * ec->I_0to) - a_to * ec->I_1to)
        + K_fto * ec->C1_to) - K_bto * ec->I_1to;
    ecR->I_2to = (((2 * b_to * ec->I_1to) / 360 - 2 * a_to * ec->I_2to) + K_f2to * ec->C2_to) - K_b2to * ec->I_2to;

    return I_to;
}

//Compute current from the IKur channel
double atrial_cell::IKur_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;
    double Ikur;

    // PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
    double fracIKurp0 = 0.437635; // Derived quantity (IKur_PKAp(baseline)/IKurtot)
    double fracIKurpISO = 0.718207; // Derived quantity (IKur_PKAp(ISO)/IKurtot)
    double a_Kur = (2.20 - 1) / (fracIKurpISO / fracIKurp0 - 1); // changed from (1.2-1)...  weijian
    double fracIKuravail = (1 - a_Kur) + a_Kur * (dependents.ph.IKur_PKAp / fracIKurp0); // +20// with 0.1 uM ISO

    double a_ur = exp((V + 58.0) / 90.0) /8 * KQ10_Kur;
    double b_ur = exp(-(V + 43.0) / 90.0)/8  * KQ10_Kur;
    double K_fur = 0.00862 * KQ10_Kur;
    double K_bur = 0.00227 * KQ10_Kur;

    double C_1ur = 1 - ec->C_2ur - ec->C_3ur - ec->C_4ur - ec->I_ur - ec->O_ur;
    ecR->C_2ur = 4 * a_ur * C_1ur - b_ur * ec->C_2ur + 2 * b_ur * ec->C_3ur - 3 * a_ur * ec->C_2ur;
    ecR->C_3ur = 3 * a_ur * ec->C_2ur - 2 * b_ur * ec->C_3ur + 3 * b_ur * ec->C_4ur - 2 * a_ur * ec->C_3ur;
    ecR->C_4ur = 2 * a_ur * ec->C_3ur - 3 * b_ur * ec->C_4ur + 4 * b_ur * ec->O_ur - a_ur * ec->C_4ur;
    ecR->O_ur = a_ur * ec->C_4ur - 4 * b_ur * ec->O_ur + K_bur * ec->I_ur - K_fur * ec->O_ur;
    ecR->I_ur = K_fur * ec->O_ur - K_bur * ec->I_ur;

    Ikur = pow(Q10Kur_conductance, Qpow) * fracIKuravail * g_Kur * ec->O_ur * (V - dependents.Ek);

    return Ikur;
}

//Compute current from the IKss channel
double atrial_cell::IKss_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;
    double I_Kss;

    double ass = 1.0 / (1.0 + (exp((-(V - 15.0) / 34.55))));
    double tau_Kss = (11.5 * (exp((-0.03 * V))) + 5.0) / KQ10_Kss;

    ecR->aKss = (ass - ec->aKss) / tau_Kss;

    I_Kss = KQ10_Kss * g_Kss * ec->aKss * (V - dependents.Ek);

    return I_Kss;
}

//Compute current from the IKr channel
double atrial_cell::IKr_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;

    double I_Kr_cur;

    double alpha0 = 0.00985942 * exp(2.56963 * V * FoRT);
    double alpha1 = 0.642537;
    double alpha2 = 0.105987 * exp(2.43874 * V * FoRT);
    double beta0 = 0.0104051 * exp(2.02278 * V * FoRT);
    double beta1 = 2.14838;
    double beta2 = 0.00330000 * exp(-0.577000 * V * FoRT);
    double alphai = 0.392662 * exp(-1.29659 * V * FoRT) + 0.221961;
    double betai = 0.696498 * exp(0.359760 * V * FoRT);
    double mu = (alphai * beta2) / betai;

    I_Kr_cur = g_Kr * ec->O_Kr * (V - dependents.Ek);

    ecR->C3_Kr = beta0 * ec->C2_Kr - alpha0 * ec->C3_Kr; //C3_Kr
    ecR->C1_Kr = (alpha1 * ec->C2_Kr + beta2 * ec->O_Kr + mu * ec->I_Kr_state) - (beta1 + 2.00000 * alpha2) * ec->C1_Kr; //C1_Kr
    ecR->C2_Kr = (beta1 * ec->C1_Kr + alpha0 * ec->C3_Kr) - (alpha1 + beta0) * ec->C2_Kr; //C2_Kr
    ecR->O_Kr = (alphai * ec->I_Kr_state + alpha2 * ec->C1_Kr) - (betai + beta2) * ec->O_Kr; //O_Kr
    ecR->I_Kr_state = (alpha2 * ec->C1_Kr + betai * ec->O_Kr) - (mu + alphai) * ec->I_Kr_state; //I_Kr

    return I_Kr_cur;
}

//Compute current from the IKACh channel
double atrial_cell::IKACh_channel()
{
  const eccs *ec = &(states.ec);
  eccs *ecR = &(states_dot.ec);

  double V = ec->V;
  double L = config.ACh;

  double I_KAch;

  double alpha_j = 73.1;   // [1/s]
  double beta_j = 120 / (1 + exp(-(V + 70) / 15));

  double alpha_k = 3.7;    // [1/s]
  double beta_k = 5.82 / (1 + exp(-(V + 70) / 15));

  double n_KACh = 1.5;
  double Km_KACh = 0.28;  // [uM]
  double K_ACh = pow(L, n_KACh) / (Km_KACh + pow(L, n_KACh));

  I_KAch = g_KACh * ec->j_KACh * ec->k_KACh * K_ACh * _Ko_/(10 + _Ko_) * (V - dependents.Ek - 6.5) / (1 + 0.5 * exp((V - dependents.Ek - 120) * FoRT / 2.5));

  ecR->j_KACh = alpha_j * (1-ec->j_KACh) - beta_j * ec->j_KACh;
  ecR->k_KACh = alpha_k * (1-ec->k_KACh) - beta_k * ec->k_KACh;

  return I_KAch;
}

//Compute current from the IK1 channel
double atrial_cell::IK1_channel()
{
    const eccs *ec = &(states.ec);

    double V = ec->V;
    double I_k1;

    I_k1 = g_K1 * (_Ko_ / (_Ko_ + 0.900))* (V - dependents.Ek - 9) / (1 + exp(0.09457 * (V - dependents.Ek - 5) + 0.42984));
    if (config._CaMKII_level_ == 1 && config.ACUTE == 0) { //chronic CaMKII effects
        I_k1 = 0.6 * I_k1;  //modified from 0.5 to 0.6  weijian
    }
    return I_k1;
}

//Compute current from the IKCa channel
double atrial_cell::IKCa_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;
    double g_KCa = 4.4e-6 * pow(V, 3) + 4e-4 * pow(V, 2) + 0.0164 * V + 0.7488;   // [mS/uF]
    double I_KCa = g_KCa * (ec->O1_KCa + ec->O2_KCa) * (V - dependents.Ek);

    double Ca_i = ec->Ca_i * 1e3;   // [uM]

    double kf1 = 120;  // [uM-1 s-1]
    double kf2 = 96;  // [uM-1 s-1]
    double kf3 = 48;  // [uM-1 s-1]
    double kb1 = 80;  // [s-1]
    double kb2 = 80;  // [s-1]
    double kb3 = 200;  // [s-1]
    double kfo1 = 160;  // [s-1]
    double kfo2 = 1200;  // [s-1]
    double kbo1 = 1000;  // [s-1]
    double kbo2 = 100;  // [s-1]

    double C1 = 1 - ec->C2_KCa - ec->C3_KCa - ec->C4_KCa - ec->O1_KCa - ec->O2_KCa;
    ecR->C2_KCa = kf1 * Ca_i * C1 + kb2 * ec->C3_KCa - (kb1 + kf2 * Ca_i) * ec->C2_KCa;
    ecR->C3_KCa = kf2 * Ca_i * ec->C2_KCa + kb3 * ec->C4_KCa + kbo1 * ec->O1_KCa - (kb2 + kf3 * Ca_i + kfo1) * ec->C3_KCa;
    ecR->C4_KCa = kf3 * Ca_i * ec->C3_KCa + kbo2 * ec->O2_KCa - (kb3 + kfo2) * ec->C4_KCa;
    ecR->O1_KCa = kfo1 * ec->C3_KCa - kbo1 * ec->O1_KCa;
    ecR->O2_KCa = kfo2 * ec->C4_KCa - kbo2 * ec->O2_KCa;

    return I_KCa;
}

//INab, I_ClCa, I_Clb, I_Cap, I_Cab
    double atrial_cell::background_currents(const int type)
    {
        const eccs *ec = &(states.ec);

        double V = ec->V;

        if (type == 1) {
            double GNaB = dependents.GNaB;
            if (config._CaMKII_level_ == 1 && config._LOOP_ == 0) GNaB = GNaB * 4; // chronic CaMKII effects
            if (config.PLMkoFlag == 1) GNaB = GNaB * 48.0 / 20.0;
            dependents.I_Nab_junc = config.INab_Block * Fjunc_na * GNaB * (V - dependents.Ena_junc);
            dependents.I_Nab_sl = config.INab_Block * Fsl_na * GNaB * (V - dependents.Ena_sl);

            double I_Nab = dependents.I_Nab_junc + dependents.I_Nab_sl;

            return I_Nab;

        } else if (type == 2) {
            double GClCa = 0.109625; // [mS/uF]
            double KdClCa = 100e-3; // [mM]

            double I_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / ec->Ca_junc) * (V - dependents.Ecl);
            double I_ClCa_sl = Fsl * GClCa / (1 + KdClCa / ec->Ca_sl) * (V -  dependents.Ecl);
            double I_ClCa = I_ClCa_junc + I_ClCa_sl;

            return I_ClCa;

        } else if (type == 3) {
            double GClB = 9e-3; // [mS/uF]
            double I_Clb = GClB * (V - dependents.Ecl);

            return I_Clb;

        } else if (type == 4) {
            double I_Cap_junc = Fjunc * _IbarpCa_ * pow(ec->Ca_junc, 1.6) / (pow(_KmPCa_, 1.6) + pow(ec->Ca_junc, 1.6));
            double I_Cap_sl = Fsl * _IbarpCa_ * pow(ec->Ca_sl, 1.6) / (pow(_KmPCa_, 1.6) + pow(ec->Ca_sl, 1.6));
            double I_Cap = I_Cap_junc + I_Cap_sl;

            dependents.I_Cap_junc = I_Cap_junc;
            dependents.I_Cap_sl = I_Cap_sl;

            return I_Cap;

        } else if (type == 5) {
            double I_Cab_junc = config.ICab_Block * Fjunc * _GCaB_ * (V - dependents.Eca_junc);
            double I_Cab_sl = config.ICab_Block * Fsl * _GCaB_ * (V - dependents.Eca_sl);
            double I_Cab = I_Cab_junc + I_Cab_sl;
            //Caffeine effect
            if (config.caffeine_on && time > config.caffeine_time) I_Cab *= 0;
            dependents.I_Cab_junc = I_Cab_junc;
            dependents.I_Cab_sl = I_Cab_sl;

            return I_Cab;

        } else if (type == 6) {
            // I_Kb, background K+
            double xkb = 1.0 / (1.0 + exp(-(V + 10) / 37.34));

            dependents.I_Kb = _GKB_ * xkb * (V - dependents.Ek);

            return dependents.I_Kb;

        } else {

            return 0;
        }
}

double atrial_cell::INCX_channel()
{
    const eccs *ec = &(states.ec);

    double V = ec->V;
    double Ca_junc = ec->Ca_junc;
    double Ca_sl = ec->Ca_sl;
    double Nai_junc = ec->Nai_junc;
    double Nai_sl = ec->Nai_sl;
    double IbarNCX = _IbarNCX_; // [uA/uF]

    if (config._CaMKII_level_ == 1 && config.ACUTE == 0) IbarNCX = 1.5 * IbarNCX; // chronic CaMKII effects

    double Kdact = _Kdact_;
    double ksat = _ksat_;
    double nu = _nu_;

    double Ka_junc = 1 / (1 + pow((Kdact / Ca_junc), 3));
    double Ka_sl = 1 / (1 + pow((Kdact / Ca_sl), 3));

    double s1_junc = exp(nu * V * FoRT) * pow(Nai_junc, 3) * _Cao_;
    double s2_junc = exp((nu - 1) * V * FoRT) * pow(_Nao_, 3) * Ca_junc;
    double s3_junc = (_KmCai_ * pow(_Nao_, 3) * (1 + pow((Nai_junc / _KmNai_), 3)) + pow(_KmNao_, 3) * Ca_junc
        + pow(_KmNai_, 3) * _Cao_ * (1 + Ca_junc / _KmCai_) + _KmCao_ * pow(Nai_junc, 3) + pow(Nai_junc, 3) * _Cao_
        + pow(_Nao_, 3) * Ca_junc) * (1 + ksat * exp((nu - 1) * V * FoRT));
    double I_NCX_junc = config.INCX_Block * Fjunc_ncx * IbarNCX * Ka_junc * (s1_junc - s2_junc) / s3_junc;

    double s1_sl = exp(nu * V * FoRT) * pow(Nai_sl, 3) * _Cao_;
    double s2_sl = exp((nu - 1) * V * FoRT) * pow(_Nao_, 3) * Ca_sl;
    double s3_sl = (_KmCai_ * pow(_Nao_, 3) * (1 + pow((Nai_sl / _KmNai_), 3)) + pow(_KmNao_, 3) * Ca_sl
        + pow(_KmNai_, 3) * _Cao_ * (1 + Ca_sl / _KmCai_)  + _KmCao_ * pow(Nai_sl, 3) + pow(Nai_sl, 3) * _Cao_
        + pow(_Nao_, 3) * Ca_sl) * (1 + ksat * exp((nu - 1) * V * FoRT));
    double I_NCX_sl = config.INCX_Block * Fsl_ncx * IbarNCX * Ka_sl * (s1_sl - s2_sl) / s3_sl;

    double I_NCX = (I_NCX_junc + I_NCX_sl);

    dependents.I_NCX_junc = I_NCX_junc;
    dependents.I_NCX_sl = I_NCX_sl;

    return I_NCX;
}

double atrial_cell::LCC_channel()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double V = ec->V;
    double Ca_junc = ec->Ca_junc;
    double Ca_sl = ec->Ca_sl;
    double Nai_junc = ec->Nai_junc;
    double Nai_sl = ec->Nai_sl;
    double LCC_CKp = dependents.ph.LCC_CKp;
    double LCCa_PKAp = dependents.ph.LCCa_PKAp;
    double LCCb_PKAp = dependents.ph.LCCb_PKAp;

    double K_Ica = 1.65; // MOUSE
    double pNa = config.ICaL_Block * K_Ica * 1.5e-8; // [cm/sec]
    double pCa = config.ICaL_Block * K_Ica * 5.4e-4; // [cm/sec] - Ca permeability
    double pK = config.ICaL_Block * K_Ica * 2.7e-7; // [cm/sec]

    double ibark = pK * (V * Frdy * FoRT) * (0.75 * ec->Ki * exp(V * FoRT) - 0.75 * _Ko_) / (exp(V * FoRT) - 1);
    double ibarna_j = pNa * (V * Frdy * FoRT) * (0.75 * Nai_junc * exp(V * FoRT) - 0.75 * _Nao_) / (exp(V * FoRT) - 1);
    double ibarna_sl = pNa * (V * Frdy * FoRT) * (0.75 * Nai_sl * exp(V * FoRT) - 0.75 * _Nao_) / (exp(V * FoRT) - 1);

    // LCC Current Fixed Parameters
    double taupo = 1.0;          // [ms] - Time constant of activation
    double TBa = 450.0;          // [ms] - Time constant
    double s1o = 0.0221;
    double k1o = 0.03;
    double kop = 2.5e-3;       // [mM]
    double cpbar = 8e-3;       // [mM]
    double tca = 78.0312;
    double ICa_scale = 5.25;
    double recoveryReduc = 3.0;

    // modified for atrial cell (weijian)
    ICa_scale = ICa_scale_ratio * ICa_scale;

    // PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
    double fracLCCbp0 = 0.250657; // Derived quantity - (LCCbp(baseline)/LCCbtot)
    double fracLCCbpISO = 0.525870; // Derived quantity - (LCCbp(ISO)/LCCbtot)
    double a_favail = (2.56 - 1) / (fracLCCbpISO / fracLCCbp0 - 1); // fracLCCbp ISO (x1.56 o.1 ISO)  //changed from (1.56-1)... (weijian)
    double favail = (1 - a_favail) + a_favail * (LCCb_PKAp / fracLCCbp0); // Test (max x2.52 100% phosph)
    ICa_scale = ICa_scale * favail;

    // CaMKII AND PKA-DEPENDENT SHIFTING OF JUNCTIONAL LCCS TO MODE 2
    double fracLCCap0 = 0.219577; // Derived
    double frac_fpkam2 = (0.15 * fracLCCap0) / (1 - fracLCCap0);
    double fpkam2 = (0.15 + frac_fpkam2) * LCCa_PKAp - frac_fpkam2; // Assumes max (100%) phos results in 15% mode 2 channels

    double fckiim2 = config.max_m2_LCC_CK * LCC_CKp; // Assumes max phos results in 10% mode 2 channels (max LCC_CKp = 1)
    // Sum up total fraction of CKII and PKA-shifted mode 2 channels
    double junc_mode2 = fckiim2 + fpkam2;

    double SSAshift = 4.4;  // modified for atrial cell (weijian)
    double SSIshift = 0;

    // Voltage- and Ca-dependent Parameters
    double poss = 1 / (1 + exp(-(V + SSAshift) / 6.14));  // modified for atrial cell (weijian)
    double Rv = 10 + 4954.0 * exp(V / 15.6);
    double PrLCC = 1 - 1 / (1 + exp(-(V + 40) / 4.0));
    double PsLCC = 1 / (1 + exp(-(V + 40 + SSIshift) / 11.32));

    double fcaj = 1 / (1 + pow((kop / Ca_junc), 3));
    double TCaj = (tca + 0.1 * (1 + pow((Ca_junc / cpbar), 2))) / (1 + pow((Ca_junc / cpbar), 2));
    double tauCaj = (Rv - TCaj) * PrLCC + TCaj;
    double tauBa = (Rv - TBa) * PrLCC + TBa;

    // Tranisition Rates (20 rates)
    double alphaLCC = poss / taupo;
    double betaLCC = (1 - poss) / taupo;
    double r1 = 0.3;                               // [1/ms] - Opening rate
    double r2 = 3.0;                                 // [1/ms] - closing rate
    double s1 = s1o * fcaj;
    double s1p = .00195;                           // [ms] - Inactivation rate
    double k1 = k1o * fcaj;
    double k1p = .00413;                           // [ms] - Inactivation rate
    double k2 = 1e-4;                              // [ms] - Inactivation rate
    double k2p = .00224;                           // [ms] - Inactivation rate
    double s2 = s1 * (k2 / k1) * (r1 / r2);
    double s2p = s1p * (k2p / k1p) * (r1 / r2);
    double k3 = exp(-(V + 40) / 3.0) / (3 * (1 + exp(-(V + 40) / 3.0)));
    double k3p = k3;
    double k5 = (1 - PsLCC) / tauCaj;
    double k6 = (fcaj * PsLCC) / tauCaj;
    double k5p = (1 - PsLCC) / tauBa;

    // Recovery terms
    k5 = k5 / recoveryReduc;
    k5p = k5p / recoveryReduc;
    double k6p = PsLCC / tauBa;
    double k4 = k3 * (alphaLCC / betaLCC) * (k1 / k2) * (k5 / k6);
    double k4p = k3p * (alphaLCC / betaLCC) * (k1p / k2p) * (k5p / k6p);

    /////////// MODE 1 junctional LCCs /////////////

    double Po_LCCj_m1 = 1.0 - ec->C2_m1j - ec->C1_m1j - ec->I1Ca_m1j - ec->I2Ca_m1j - ec->I1Ba_m1j - ec->I2Ba_m1j;

    double LCC_m1[] = {Po_LCCj_m1, ec->C1_m1j, ec->C2_m1j, ec->I1Ca_m1j, ec->I2Ca_m1j, ec->I1Ba_m1j, ec->I2Ba_m1j,
        alphaLCC, betaLCC, r1, r2, s1, s2, k1, k2, k3, k4, k5, k6, s1p, s2p, k1p, k2p, k3p, k4p, k5p, k6p};

    LCC_rates(LCC_m1, &(ecR->C2_m1j));

    double ibarca_jm1 = (4 * pCa * V * Frdy * FoRT) * (0.001 * exp(2 * V * FoRT) - 0.341 * _Cao_) / (exp(2 * V * FoRT) - 1);
    double I_Ca_junc_m1 = (Fjunc_CaL * ibarca_jm1 * Po_LCCj_m1) * ICa_scale;

    ////////////

    /////////// MODE 2 junctional LCCs /////////////

    // Re-define parameters as mode 2 specific parameters
    double r2m2 = 3.0 / 8.0; // [1/ms] - closing rate,  changed from rabbit (3/10) - MOUSE
    double s2m2 = s1 * (k2 / k1) * (r1 / r2m2);
    double s2pm2 = s1p * (k2p / k1p) * (r1 / r2m2);

    // State transitions for MODE 2 junctional LCCs
    double Po_LCCj_m2 = 1.0 - ec->C2_m2j - ec->C1_m2j - ec->I1Ca_m2j - ec->I2Ca_m2j - ec->I1Ba_m2j - ec->I2Ba_m2j; // O_m2j

    double LCC_m2[] = {Po_LCCj_m2, ec->C1_m2j, ec->C2_m2j, ec->I1Ca_m2j, ec->I2Ca_m2j, ec->I1Ba_m2j, ec->I2Ba_m2j,
        alphaLCC, betaLCC, r1, r2m2, s1, s2m2, k1, k2, k3, k4, k5, k6, s1p, s2pm2, k1p, k2p, k3p, k4p, k5p, k6p};

    LCC_rates(LCC_m2, &(ecR->C2_m2j));

    double ibarca_jm2 = (4 * pCa * V * Frdy * FoRT) * (.001 * exp(2 * V * FoRT) - 0.341 * _Cao_) / (exp(2 * V * FoRT) - 1);
    double I_Ca_junc_m2 = (Fjunc_CaL * ibarca_jm2 * (Po_LCCj_m2)) * ICa_scale;

    ///////////

    // Total junctional ICa
    double I_Ca_junc = (1 - junc_mode2) * I_Ca_junc_m1 + junc_mode2 * I_Ca_junc_m2;

    // SUB-SARCOLEMMAL LCCs
    // Re-assign necessary params to be Casl sensitive
    double fcasl = 1 / (1 + pow((kop / Ca_sl), 3));    // Depends on sl Ca
    double TCasl = (tca + 0.1 * pow((1 + (Ca_sl / cpbar)), 2)) / (1 + pow((Ca_sl / cpbar), 2));
    double tauCasl = (Rv - TCasl) * PrLCC + TCasl;

    // Re-assign necessary rates to be Casl sensitive
    double s1sl = s1o * fcasl;
    double k1sl = k1o * fcasl;
    double s2sl = s1sl * (k2 / k1sl) * (r1 / r2);
    double s2psl = s1p * (k2p / k1p) * (r1 / r2);
    double k5sl = (1 - PsLCC) / tauCasl;
    k5sl = k5sl / recoveryReduc;  // Reduced for recovery
    double k6sl = (fcasl * PsLCC) / tauCasl;
    double k4sl = k3 * (alphaLCC / betaLCC) * (k1sl / k2) * (k5sl / k6sl);
    double k4psl = k3p * (alphaLCC / betaLCC) * (k1p / k2p) * (k5p / k6p);

    /////////// MODE 1 sarcolemmal LCCs /////////////

    // State transitions for 'mode 1' sarcolemmal LCCs
    double Po_LCCsl_m1 = 1 - ec->C2_m1sl - ec->C1_m1sl - ec->I1Ca_m1sl - ec->I2Ca_m1sl - ec->I1Ba_m1sl - ec->I2Ba_m1sl; // O_m1sl
    double LCC_slm1[] = {Po_LCCsl_m1, ec->C1_m1sl, ec->C2_m1sl, ec->I1Ca_m1sl, ec->I2Ca_m1sl, ec->I1Ba_m1sl, ec->I2Ba_m1sl,
        alphaLCC, betaLCC, r1, r2, s1sl, s2sl, k1sl, k2, k3, k4sl, k5sl, k6sl, s1p, s2psl, k1p, k2p, k3p, k4psl, k5p, k6p};

    LCC_rates(LCC_slm1, &(ecR->C2_m1sl));

    double ibarca_slm1 = (4 * pCa * V * Frdy * FoRT) * (.001 * exp(2 * V * FoRT) - 0.341 * _Cao_) / (exp(2 * V * FoRT) - 1);
    double I_Casl_m1 = (Fsl_CaL * ibarca_slm1 * Po_LCCsl_m1) * ICa_scale;

    ///////////////

    /////////// MODE 2 sarcolemmal LCCs /////////////
    double r2slm2 = r2m2;
    double s2slm2 = s1sl * (k2 / k1sl) * (r1 / r2slm2);
    double s2pslm2 = s1p * (k2p / k1p) * (r1 / r2slm2);

    // State transitions for mode 2 sarcolemmal LCCs
    double Po_LCCsl_m2 = 1 - ec->C2_m2sl - ec->C1_m2sl - ec->I1Ca_m2sl- ec->I2Ca_m2sl - ec->I1Ba_m2sl - ec->I2Ba_m2sl; //O_m2sl

    double LCC_slm2[27] = {Po_LCCsl_m2, ec->C1_m2sl, ec->C2_m2sl, ec->I1Ca_m2sl, ec->I2Ca_m2sl, ec->I1Ba_m2sl, ec->I2Ba_m2sl,
        alphaLCC, betaLCC, r1, r2slm2, s1sl, s2slm2, k1sl, k2, k3, k4sl, k5sl, k6sl, s1p, s2pslm2, k1p, k2p, k3p, k4psl, k5p, k6p};

    LCC_rates(LCC_slm2, &(ecR->C2_m2sl));

    double ibarca_slm2 = (4 * pCa * V * Frdy * FoRT) * (.001 * exp(2 * V * FoRT) - 0.341 * _Cao_) / (exp(2 * V * FoRT) - 1);
    double I_Casl_m2 = (Fsl_CaL * ibarca_slm2 * Po_LCCsl_m2) * ICa_scale;

    ///////////////

    // Sum mode 1 and mode 2 sl channels for total sl current
    double fckiim2_sl = 0; // Set to zero since SL LCCp by CaMKII is negligible
    double sl_mode2 = fckiim2_sl + fpkam2;
    double I_Ca_sl = (1 - sl_mode2) * I_Casl_m1 + sl_mode2 * I_Casl_m2;

    // Total currents through LCC
    double I_Ca = I_Ca_junc + I_Ca_sl;

    double I_CaK_junc = ibark * Fjunc_CaL * ((1 - junc_mode2) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2) * ICa_scale;
    double I_CaK_sl = ibark * Fsl_CaL * ((1 - sl_mode2) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * ICa_scale;

    double I_CaNa_junc = (Fjunc_CaL * ibarna_j * ((1 - junc_mode2) * Po_LCCj_m1 + junc_mode2 * Po_LCCj_m2)) * ICa_scale;
    double I_CaNa_sl = Fsl_CaL * ibarna_sl * ((1 - sl_mode2) * Po_LCCsl_m1 + sl_mode2 * Po_LCCsl_m2) * ICa_scale;

    double I_Ca_tot_junc = I_Ca_junc + dependents.I_Cab_junc + dependents.I_Cap_junc - 2 * dependents.I_NCX_junc; // [uA/uF]
    double I_Ca_tot_sl = I_Ca_sl + dependents.I_Cab_sl + dependents.I_Cap_sl - 2 * dependents.I_NCX_sl;   // [uA/uF]

    dependents.I_Ca_junc = I_Ca_junc;
    dependents.I_Ca_sl = I_Ca_sl;
    dependents.I_CaNa_junc = I_CaNa_junc;
    dependents.I_CaNa_sl = I_CaNa_sl;
    dependents.I_CaK_junc = I_CaK_junc;
    dependents.I_CaK_sl = I_CaK_sl;
    dependents.I_Ca_tot_junc = I_Ca_tot_junc;
    dependents.I_Ca_tot_sl = I_Ca_tot_sl;

    return I_Ca;
}

void atrial_cell::LCC_rates(const double *p, double *ptr)
{
    double O = p[0];
    double C1 = p[1];
    double C2 = p[2];
    double I1Ca = p[3];
    double I2Ca = p[4];
    double I1Ba = p[5];
    double I2Ba = p[6];
    double alpha = p[7];
    double beta = p[8];
    double r1 = p[9];
    double r2 = p[10];
    double s1 = p[11];
    double s2 = p[12];
    double k1 = p[13];
    double k2 = p[14];
    double k3 = p[15];
    double k4 = p[16];
    double k5 = p[17];
    double k6 = p[18];
    double s1p = p[19];
    double s2p = p[20];
    double k1p = p[21];
    double k2p = p[22];
    double k3p = p[23];
    double k4p = p[24];
    double k5p = p[25];
    double k6p = p[26];

    *ptr = beta * C1 + k5 * I2Ca + k5p * I2Ba - (k6 + k6p + alpha) * C2;
    *(++ptr) = alpha * C2 + k2 * I1Ca + k2p * I1Ba + r2 * O - (r1 + beta + k1 + k1p) * C1;
    *(++ptr) = k1 * C1 + k4 * I2Ca + s1 * O - (k2 + k3 + s2) * I1Ca;
    *(++ptr) = k3 * I1Ca + k6 * C2 - (k4 + k5) * I2Ca;
    *(++ptr) = k1p * C1 + k4p * I2Ba + s1p * O - (k2p + k3p + s2p) * I1Ba;
    *(++ptr) = k3p * I1Ba + k6p * C2 - (k5p + k4p) * I2Ba;
}

void atrial_cell::Na_concentrations()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    double Nai_junc = ec->Nai_junc;
    double Nai_sl = ec->Nai_sl;

    // Na buffers
    ecR->NaBj = _kon_na_ * Nai_junc * (_Bmax_Naj_ - ec->NaBj) - _koff_na_ * ec->NaBj;  // NaBj      [mM/ms]
    ecR->NaBsl = _kon_na_ * Nai_sl * (_Bmax_Nasl_ - ec->NaBsl) - _koff_na_ * ec->NaBsl; // NaBsl     [mM/ms]

    // Na Concentrations
    double I_Na_tot_junc = dependents.I_Na_junc + dependents.I_Nab_junc + 3 * dependents.I_NCX_junc
        + 3 * dependents.I_Nak_junc + dependents.I_CaNa_junc; // [uA/uF]
    double I_Na_tot_sl = dependents.I_Na_sl + dependents.I_Nab_sl + 3 * dependents.I_NCX_sl
        + 3 * dependents.I_Nak_sl + dependents.I_CaNa_sl;   //[uA/uF]

    ecR->Nai_junc = -I_Na_tot_junc * Cmem / (Vjunc * Frdy) + J_na_juncsl / Vjunc * (Nai_sl - Nai_junc) - ecR->NaBj; //Nai_junc
    ecR->Nai_sl = -I_Na_tot_sl * Cmem / (Vsl * Frdy) + J_na_juncsl / Vsl * (Nai_junc - Nai_sl) - ecR->NaBsl
        + J_na_slmyo / Vsl * (ec->Nai_i - Nai_sl); //Nai_sl

    if (config.NaClampFlag == 1) ecR->Nai_i = 0; // Na clamp
    else {
        ecR->Nai_i = J_na_slmyo / Vmyo * (Nai_sl - ec->Nai_i); // [mM/msec], Nai_i
    }
}

void atrial_cell::classic_Ca_handling()
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);
    
    double Cai = ec->Ca_i;
    double Ca_junc = ec->Ca_junc;
    double Ca_sl = ec->Ca_sl;
    double Ca_sr = ec->Ca_sr;
    double Csqn = ec->Csqn;

    // PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
    double fracTnIpo = 0.062698; // Derived quantity (TnI_PKAp(baseline)/TnItot)
    double fPKA_TnI = (1.61 - 0.61 * (1 - dependents.ph.TnI_PKAp) / (1 - fracTnIpo)); // Max effect +61//    //fPKA_TnI = (1.45-0.45*(1-TnI_PKAp)/(1-fracTnIpo)); // Max effect +45//
    double koff_tncl = _koff_tncl_ * fPKA_TnI;

    // IP3r Update
    if (config.IP3R_effect_on) {
        double pIP3r[] = {ec->IP3r_Junc_O, ec->IP3r_Junc_RI, ec->IP3r_Junc_R, ec->IP3r_Junc_RC, ec->IP3r_Junc_RC2,
            ec->IP3r_Junc_RC3, Cai * 1e3, ec->IP3_state};
        IP3R_rates(pIP3r, &(ecR->IP3r_Junc_O));
    }

    double pRyR[] = {ec->Ca_sr, ec->Ca_junc, ec->R_SR, ec->O_SR, ec->I_SR, ec->IP3r_Junc_O};
    double J_SRCarel = RyR_rel(pRyR, &(ecR->R_SR));

    double pRyr_leak[] = {ec->Ca_sr, ec->Ca_junc};
    double J_SRleak = RyR_leak(pRyr_leak);

    double pSerca[] = {ec->Ca_sr, ec->Ca_i};
    double J_serca = SERCA(pSerca);

    // Junctional and SL Ca Buffers
    ecR->SLLj = _kon_sll_ * Ca_junc * (_Bmax_SLlowj_ - ec->SLLj) - _koff_sll_ * ec->SLLj;      // [mM/ms]
    ecR->SLHj = _kon_slh_ * Ca_junc * (_Bmax_SLhighj_ - ec->SLHj) - _koff_slh_ * ec->SLHj;      // [mM/ms]
    double J_CaB_junction = ecR->SLLj + ecR->SLHj;

    ecR->SLLsl = _kon_sll_ * Ca_sl * (_Bmax_SLlowsl_ - ec->SLLsl) - _koff_sll_ * ec->SLLsl;     // [mM/ms]
    ecR->SLHsl = _kon_slh_ * Ca_sl * (_Bmax_SLhighsl_ - ec->SLHsl) - _koff_slh_ * ec->SLHsl;     // [mM/ms]
    double J_CaB_sl = ecR->SLLsl + ecR->SLHsl;

    // Csqn      [mM/ms]
    ecR->Csqn = (_kon_csqn_ * Ca_sr * (_Bmax_Csqn_ - Csqn) - _koff_csqn_ * Csqn);

    // Ca_sr     [mM/ms] // Ratio 3 leak current
    ecR->Ca_sr = J_serca * Vmyo / Vjsr - J_SRleak * Vmyo / Vjsr - ecR->Csqn - J_SRCarel;

    double K_ip3_diffusion = 1.0;
    // if (if_ip3r_effect) {
    //   double O_IP3r = ec->IP3r_Junc_O;
    //   K_ip3_diffusion = K_ip3_diffusion * O_IP3r * 12;
    // }

    // Ca_sl
    ecR->Ca_sl += -dependents.I_Ca_tot_sl * Cmem / (Vsl * 2 * Frdy) + J_ca_juncsl / Vsl * (Ca_junc - Ca_sl) 
        - K_ip3_diffusion * J_ca_slmyo / Vsl * (Ca_sl - Cai) - J_CaB_sl;

    // Ca_junc
    ecR->Ca_junc += -dependents.I_Ca_tot_junc * Cmem / (Vjunc * 2 * Frdy) - J_CaB_junction - J_ca_juncsl / Vjunc * (Ca_junc - Ca_sl)
        + J_SRCarel * Vjsr / Vjunc + J_SRleak * Vmyo / Vjunc;

    // Cytosolic Ca Buffers
    ecR->TnCL = _kon_tncl_ * Cai * (_Bmax_TnClow_ - ec->TnCL) - koff_tncl * ec->TnCL;            // TnCL      [mM/ms]
    ecR->TnCHc = _kon_tnchca_ * Cai * (_Bmax_TnChigh_ - ec->TnCHc - ec->TnCHm) - _koff_tnchca_ * ec->TnCHc; // TnCHc     [mM/ms]
    ecR->TnCHm = _kon_tnchmg_ * _Mgi_ * (_Bmax_TnChigh_ - ec->TnCHc - ec->TnCHm) - _koff_tnchmg_ * ec->TnCHm;   // TnCHm     [mM/ms]

    ecR->Myosin_ca = _kon_myoca_ * Cai * (_Bmax_myosin_ - ec->Myosin_ca - ec->Myosin_mg) - _koff_myoca_ * ec->Myosin_ca; // Myosin_ca [mM/ms]
    ecR->Myosin_mg = _kon_myomg_ * _Mgi_ * (_Bmax_myosin_ - ec->Myosin_ca - ec->Myosin_mg) - _koff_myomg_ * ec->Myosin_mg; // Myosin_mg [mM/ms]
    ecR->SRB = _kon_sr_ * Cai * (_Bmax_SR_ - ec->SRB) - _koff_sr_ * ec->SRB; // SRB       [mM/ms]

    double J_CaB_cytosol = ecR->TnCL + ecR->TnCHc + ecR->Myosin_ca + ecR->SRB;

    // Cai
    ecR->Ca_i += -J_serca - J_CaB_cytosol + K_ip3_diffusion * J_ca_slmyo / Vmyo * (Ca_sl - Cai);

    dependents.Ca_sr = Ca_sr;
    dependents.Csqn = Csqn;
    dependents.J_CaB_junction = J_CaB_junction;
    dependents.J_CaB_sl = J_CaB_sl;
    dependents.J_CaB_cytosol = J_CaB_cytosol;
    dependents.O_SR = ec->O_SR;
    dependents.O_IP3r = ec->IP3r_Junc_O;
}

void atrial_cell::IP3R_rates(const double *p, double *ptr)
{
    double O = p[0];
    double RI = p[1];
    double r = p[2];
    double RC = p[3];
    double RC2 = p[4];
    double RC3 = p[5];
    double RC4 = 1 - O - RI - r - RC - RC2 - RC3;

    double Cai = p[6]; //[uM]
    double IP3 = p[7]; //[uM]


    double kf1 = 27.732563188;
    double kf2 = 0.9530238909;
    double kf3 = 0.0093309470;
    double kf4 = 0.0014238430;
    double kf5 = 0.2099957248;
    double kf6 = 11.564288084;
    double kr1 = 1.5785758365;
    double kr2 = 0.1047000479;
    double kr3 = 0.0005794316;
    double kr4 = 1.1971936265;
    double kr5 = 2.1134860307;
    double kr6 = 7.0688294419;

    *ptr = RI * Cai * kf1 - O * kr1; //O
    *(++ptr) = (O * kr1 + r * IP3 * kf2) - RI * (Cai * kf1 + kr2);    //RI
    *(++ptr) = (RI * kr2 + RC * kr3) - r * (IP3 * kf2 + Cai * kf3);        //r
    *(++ptr) = (r * Cai * kf3 + RC2 * kr4) - RC * (kr3 + Cai * kf4);    //RC
    *(++ptr) = (RC * Cai * kf4 + RC3 * kr5) - RC2 * (kr4 + Cai * kf5);    //RC2
    *(++ptr) = (RC2 * Cai * kf5 + RC4 * kr6) - RC3 * (kr5 + Cai * kf6); //RC3
}

double atrial_cell::RyR_rel(const double *p, double *ptr)
{
    const double Ca_sr = p[0];
    const double Ca_junc = p[1];
    const double r = p[2];
    const double O = p[3];
    const double I = p[4];
    const double O_IP3r = p[5];

    // PKA phosphoregulation
    double frac_RyRo = 0.204276; // Derived (RyR_PKAp(basal)/_RyRtot_)
    double a_RyR = (2 - 1) / (1.0 / frac_RyRo - 1); // Max effect: fPKA_RyR=2
    double fPKA_RyR = 1 - a_RyR + a_RyR * (dependents.ph.RyR_PKAp / frac_RyRo); // 1 with NO ISO

    // CaMKII phosphoregulation
    double fCKII_RyR = (10 * dependents.ph.RyR_CKp - 1); // 1 at basal condition - MOUSE
    double fCKII_ec50SR = 1.16 - 4.0 / 5.0 * dependents.ph.RyR_CKp;
    double ec50SR = fCKII_ec50SR * _ec50SR_;

    // RyR Parameters
    double MaxSR = 15.0;
    double MinSR = 1.0;
    double kCaSR = MaxSR - (MaxSR - MinSR) / (1 + pow((ec50SR / Ca_sr), 2.5));
    double kiSRCa = _kiCa_ * kCaSR;
    double koSRCa = _koCa_ / kCaSR;

    // Caffeine Inducing
    if (config.caffeine_on && time > config.caffeine_time) {
        koSRCa *= 7.5;
    }

    // Overall phosphoregulation on koSRCa
    koSRCa = (fCKII_RyR + fPKA_RyR - 1) * koSRCa;

    // ODEs for RyR states and SR release through open RyRs
    const double RI = 1 - r - O - I;

    *ptr = (_kim_ * RI - kiSRCa * Ca_junc * r) - (koSRCa * pow(Ca_junc, 2) * r - _kom_ * O); // R
    *(++ptr) = (koSRCa * pow(Ca_junc, 2) * r - _kom_ * O) - (kiSRCa * Ca_junc * O - _kim_ * I);    // O
    *(++ptr) = (kiSRCa * Ca_junc * O - _kim_ * I) - (_kom_ * I - koSRCa * pow(Ca_junc, 2) * RI);   // I

    double J_rel = _ks_ * O * (Ca_sr - Ca_junc);  // [mmol/L SR/ ms]

    dependents.J_SRCarel = J_rel;

    return J_rel;
}

double atrial_cell::RyR_leak(const double *p)
{
    const eccs *ec = &(states.ec);

    double Ca_sr = p[0];
    double Ca_junc = p[1];

    double kleak = 2 * 5.348e-6; // [1/ms] changed from rabbit (5.348e-6)

    // Passive RyR leak - includes CaMKII regulation of leak flux
    kleak = (0.5 + 2.5 * dependents.ph.RyR_CKp) * kleak; // MOUSE (reduced CaMKII effect on leak)
    double J_SRleak = config.JSRleak_Block * kleak * (Ca_sr - Ca_junc); // [mmol/L cyt/ms]

    ///////// IP3R //////////
    // if (if_ip3r_effect) //IP3 Effects
    // {
    //   kleak /= 2;
    //   double Kleak_ip3 = kleak * 12 * ec->IP3r_Junc_O;
    //
    //   double J_IP3r = Kleak_ip3 * (Ca_sr - Ca_junc);
    //   J_SRleak += J_IP3r;
    //
    //   sh->J_IP3r = J_IP3r;
    // }

    dependents.J_SRleak = J_SRleak;

    return J_SRleak;
}

double atrial_cell::SERCA(const double *p)
{
    double Ca_sr = p[0];
    double Cai = p[1];

    // CaMKII-dependent phosphoregulation
    // double fCKII_PLB = (1-0.5*ph->PLB_CKp); // Max effect: fCKII_PLB=0.5
    double fCKII_PLB = (1 - 0.5 * dependents.ph.PLB_CKp);

    // PKA-dependent phosphoregulation
    double fracPKA_PLBo = 1 - 0.079755; // Derived quantity - (1 - (PLBp(baseline)/_PLBtot_))
    double fPKA_PLB = (dependents.ph.PLB_PKAn / fracPKA_PLBo) * 1.0 / 4.0 + 3.0 / 4.0; // Max effect: fPKA_PLB=0.25
    // double fPKA_PLB = (ph->PLB_PKAn/fracPKA_PLBo)*(100-55.31)/100.0 + 55.31/100.0; // Max effect: fPKA_PLB=0.45

    double Kmf = _Kmf_;
    double Kmr = _Kmr_;

    double Kserca = 1;
    if (config.serca_amp) {
        // CaMKII effects on SERCA, Vmax=697 and 1486 μM/s at 0.2 and 4 Hz, respectively (Picht 2006), 2.132 times
        Kserca *= 1 / (pow(0.0167 / (dependents.ph.PLB_CKp), 12) + 1) + 1;
    }

    //IP3 Effects
    // if (if_ip3r_effect)
    // {
    //     Kmf = Kmf * (1.0+0.4*K_IP3_up);
    // }

    // Select smaller value (resulting in max reduction of Kmf)
    double fKmf = (fCKII_PLB <= fPKA_PLB) ? fCKII_PLB : fPKA_PLB;

    Kmf = fKmf * Kmf;

    double J_serca = config.Jserca_Block * Kserca * _Vmax_SERCA_ * (pow((Cai / Kmf), _hillSRCaP_) - pow((Ca_sr / Kmr), _hillSRCaP_))
            / (1 + pow((Cai / Kmf), _hillSRCaP_) + pow((Ca_sr / Kmr), _hillSRCaP_)); // [mM/ms]

    if (config.caffeine_on && time > config.caffeine_time) {
        J_serca *= 0;
    }

    dependents.J_serca = J_serca;

    return J_serca;
}