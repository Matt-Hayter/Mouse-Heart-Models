#ifndef ATRIAL_CELL_HPP
#define ATRIAL_CELL_HPP

#define _USE_MATH_DEFINES //For M_P

//Atrial cell class

#include "base_cell.hpp"
#include "simulation_configuration.hpp"
#include "atrial_cell_configuration.hpp"
#include "SAN_cell.hpp"

#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>

class atrial_cell : public base_cell
{
    protected:

        //Model state variables:

        //eccs
        struct eccs 
        {
            double V;                   //1 [mV]                                                                               
            double Nai_junc;            //2 [mM]
            double Nai_sl;              //3 [mM]
            double Nai_i;               //4 [mM]
            double Ki;                  //5 [mM]
            double Ca_junc;             //6 [mM]
            double Ca_sl;               //7 [mM]
            double Ca_i;                //8 [mM]
            double Ca_sr;               //9 [mM]
            double Cl_i;                //10 [mM] !!!
            double Csqn;                //11
            double m_na;                //12 INa
            double h_na;                //13
            double j_na;                //14
            double h_l;                 //15 INaL
            double R_SR;                //16 SR
            double O_SR;                //17
            double I_SR;                //18
            double NaBj;                //19 Na Buffer
            double NaBsl;               //20
            double TnCL;                //21 Ca Buffer
            double TnCHc;               //22
            double TnCHm;               //23
            double Myosin_ca;           //24
            double Myosin_mg;           //25
            double SRB;                 //26
            double SLLj;                //27
            double SLLsl;               //28
            double SLHj;                //29
            double SLHsl;               //30
            double C2_m1j;              //31 LCC-Markov
            double C1_m1j;              //32
            double I1Ca_m1j;            //33
            double I2Ca_m1j;            //34
            double I1Ba_m1j;            //35
            double I2Ba_m1j;            //36
            double C2_m2j;              //37
            double C1_m2j;              //38
            double I1Ca_m2j;            //39
            double I2Ca_m2j;            //40
            double I1Ba_m2j;            //41
            double I2Ba_m2j;            //42
            double C2_m1sl;             //43
            double C1_m1sl;             //44
            double I1Ca_m1sl;           //45
            double I2Ca_m1sl;           //46
            double I1Ba_m1sl;           //47
            double I2Ba_m1sl;           //48
            double C2_m2sl;             //49
            double C1_m2sl;             //50
            double I1Ca_m2sl;           //51
            double I2Ca_m2sl;           //52
            double I1Ba_m2sl;           //53
            double I2Ba_m2sl;           //54
            //Added for Atria
            double C_2ur;               //55 IKur-Markov
            double C_3ur;               //56
            double C_4ur;               //57
            double O_ur;                //58
            double I_ur;                //59
            double aKss;                //60 IKss-Atria-HH
            double C3_Kr;               //61 IKr-Markov
            double C1_Kr;               //62
            double C2_Kr;               //63
            double O_Kr;                //64
            double I_Kr_state;          //65
            double C2_to;               //66 Ito-Markov
            double C1_to;               //67
            double O_to;                //68
            double I_0to;               //69
            double I_1to;               //70
            double I_2to;               //71
            double IP3r_Junc_O;         //72 IP3r
            double IP3r_Junc_RI;        //73
            double IP3r_Junc_R;         //74
            double IP3r_Junc_RC;        //75
            double IP3r_Junc_RC2;       //76
            double IP3r_Junc_RC3;       //77
            double IP3_state;           //78
            double Ca_gap;              //79 Gap Area
            double Ca_gapsr;            //80
            double Csqn_gap;            //81
            double Ca_ss_gap;           //82
            double j_KACh;              //83
            double k_KACh;              //84
            double C2_KCa;              //85
            double C3_KCa;              //86
            double C4_KCa;              //87
            double O1_KCa;              //88
            double O2_KCa;              //89
        };
        //cams: 90-104 (cj), 105-119 (cs), 120-134 (cc)
        struct cams
        {
            double CaM;                 //   CaM- Ca free CaM
            double Ca2CaM;              //   2 Ca bound to C terminal sites
            double Ca4CaM;              //   4 Ca bound
            double CaMB;                //  
            double Ca2CaMB;             //  
            double Ca4CaMB;             //  
            double Pb2;                 //   probability of a Ca2CaM bound CaMKII subunit
            double Pb;                  //   probability of a Ca4CaM bound CaMKII subunit
            double Pt;                  //   probability of a Ca4CaM bound autophosphorylated CaMKII subunit
            double Pt2;                 //   probability of a Ca2CaM bound autophosphorylated CaMKII subunit
            double Pa;                  //   probability of an autonomous autophosphorylated CaMKII subunit
            double Ca4CaN;              //  
            double CaMCa4CaN;           //   active calcineurin
            double Ca2CaMCa4CaN;        //  
            double Ca4CaMCa4CaN;        //  
        };
        //ckiis
        struct ckiis 
        {
            double LCC_CKjuncp;         //135
            double RyR2815p;            //136
            double PLBT17p;             //137
            double LCC_CKslp;           //138
        };
        //bars
        struct bars 
        {
            double LR;                  //139
            double LRG;                 //140
            double RG;                  //141
            double b1AR_S464;           //142
            double b1AR_S301;           //143
            double GsaGTPtot;           //144
            double GsaGDP;              //145
            double Gsby;                //146
            double AC_GsaGTP;           //147
            double PDEp;                //148
            double cAMPtot;             //149
            double RC_I;                //150
            double RCcAMP_I;            //151
            double RCcAMPcAMP_I;        //152
            double RcAMPcAMP_I;         //153
            double PKACI;               //154
            double PKACI_PKI;           //155
            double RC_II;               //156
            double RCcAMP_II;           //157
            double RCcAMPcAMP_II;       //158
            double RcAMPcAMP_II;        //159
            double PKACII;              //160
            double PKACII_PKI;          //161
            double I1p_PP1;             //162
            double I1ptot;              //163
            double PLBp;                //164
            double PLMp;                //165
            double LCCap;               //166
            double LCCbp;               //167
            double RyRp;                //168
            double TnIp;                //169
            double KURp;                //170
        };
        //ip3s
        struct ip3s
        {
            double Gd;                  //171 GaGDP density, [um^2-1]
            double Gt;                  //172 GaGTP density, [um^2-1]
            double r;                   //173 Noncoupled receptor density, [um^2-1]
            double Rl;                  //174 Ligand-bound receptor density, [um^2-1]
            double Rg;                  //175 Precoupled receptor density, [um^2-1]
            double Rlg;                 //176 Active receptor density, [um^2-1]
            double Rlgp;                //177 Phosphorylated receptor density, [um^2-1]
            double Pc;                  //178 PLCb-Ca density, [um^2-1]
            double Pcg;                 //179 PLCB-Ca-GaGTP density, [um^2-1]
            double P;                   //180 PLCB density, [um^2]
            double Pg;                  //181 PLCB-GaGTP density, [um^2-1]
        };

        //Holds each state type from within the state vector
        struct state_variables
        {
            eccs  ec;
            cams  cj;
            cams  cs;
            cams  cc;
            ckiis ck;
            bars  ba;
            ip3s  ip3;
            const int number_of_states{181}; //Total number of cell state variables
        };

        //Phosph
        struct phosph
        {
            //CaMKII phosphorylated targets
            double LCC_CKp; //Fractional CaMKII-dependent LCC junc phos
            double RyR_CKp; //Fractional CaMKII-dependent RyR phos
            double PLB_CKp; //Fractional CaMKII-dependent PLB phos
            //PKA phosphorylated targets
            double LCCa_PKAp;
            double LCCb_PKAp;
            double PLB_PKAn;
            double RyR_PKAp;
            double TnI_PKAp;
            double IKur_PKAp;
            double PLM_PKAp;
        };

        //Variables to be outputted
        struct outputs
        {
            //Main Currents
            double I_Na; //[pA/pF] For all
            double I_NaK;
            double I_to;
            double I_Kur;
            double I_Ks;
            double I_Kss;
            double I_Kr;
            double I_KACh;
            double I_K1;
            double I_KCa;
            double I_Kb;
            double I_NCX;
            double I_Nab;
            double I_ClCa;
            double I_Clb;
            double I_Cap;
            double I_Cab;
            double I_CaL;
            double I_j; //Junction current

            //Other currents (grouped by type)
            double I_Na_junc, I_Na_sl;
            double I_Nab_junc, I_Nab_sl;
            double I_Nak_junc, I_Nak_sl;
            double I_Ca_junc, I_Ca_sl, I_CaNa_junc, I_CaNa_sl, I_CaK_junc, I_CaK_sl;
            double I_NCX_junc, I_NCX_sl;
            double I_Cap_junc, I_Cap_sl;
            double I_Cab_junc, I_Cab_sl;
            //double I_CaT;
            double I_Ca_tot_junc, I_Ca_tot_sl;
            double I_Na_fast, I_Na_slow;
            //double I_f, I_f_hcn1, I_f_hcn2, I_f_hcn4;

            //Nerst Potentials
            double Ena_junc;
            double Ena_sl;
            double Ek, Eks;
            double Eca_junc;
            double Eca_sl;
            double Ecl;

            //Conductance
            double GNaB;
            
            //Ca2+ Handling
            double J_SRCarel;
            double J_serca;
            double J_SRleak;
            double J_CaB_junction;
            double J_CaB_sl;
            double J_CaB_cytosol;
            double Ca_sr;
            double Csqn;
            double O_SR;
            double O_IP3r;
            double J_IP3r;

            phosph ph; 

            double camkii_junc_act, camkii_cyt_act, camkii_sl_act;

            double dVdt;
        };

        //Declare cell states, currents and rates of change of states
        const simulation_configuration sim_config;
        const atrial_configuration atrial_config; //Contains user's configuration
        state_variables states; //For holding each of the state variables
        state_variables states_dot; //For holding the rates of change of each state variables
        outputs dependents; //Dependent variables to be outputted

        ////////Declare and initialise parameters and constants////////

        //Define constants
        const double _R   = 8314.0; //[J/kmol*K]
        const double Frdy = 96485.0; //[C/mol]
        const double FoRT = Frdy / _R / sim_config.temp;

        //Cell geometry
        const double junctionLength = 15e-3; //junc length [µm]
        const double junctionRadius = 160e-3; //junc radius [µm]
        const double cellL = 90.0; //cell length [µm]
        const double cellR = 6.5; //cell radius [µm]
        const double Cmem = 50.0e-12; //[F] //#define Cmem (Acell*1e-14) //[F]
        const double Vcell = M_PI * pow(cellR, 2) * cellL * 1e-15; //[L]
        const double Acell = Vcell * 0.694 * 1e15; //[µm^2]
        const double Vmyo = 0.79 * Vcell;
        const double Vjsr = 0.0176 * Vcell;
        const double Vsl = 0.017 * Vcell;
        const double Vjunc = 0.0001997 * Vcell;

        //Fractional currents
        const double Fjunc = 0.1875;
        const double Fsl = 1 - Fjunc;
        const double Fjunc_nak = 0.2268;
        const double Fsl_nak = 1 - Fjunc_nak;
        const double Fjunc_ncx = Fjunc;
        const double Fsl_ncx = 1 - Fjunc_ncx;
        const double Fjunc_CaL = 0.9;
        const double Fsl_CaL = 1 - Fjunc_CaL;
        const double Fjunc_na = Fjunc;
        const double Fsl_na = 1 - Fjunc_na;
        const double SAsl = Fsl * Acell; //[µm^2]
        const double Njunc = (Fjunc * Acell) / (M_PI * pow(junctionRadius, 2)); //[-]
        const double SAjunc = Njunc * M_PI * 2 * junctionLength * junctionRadius; //[µm^2]

        //Diffusion
        const double distSLcyto = 0.45; //distance of SL to cytosol [µm]
        const double dist_JuncSL = 0.3; //distance of junc to SL [µm]
        const double D_ca_JuncSL = 1.64e-6; //Dca junc to SL [cm^2/sec]
        const double DcaSLcyto = 1.22e-6; //Dca SL to cyto [cm^2/sec] !!!! shanzhuo changes it back
        const double D_na_JuncSL = 1.09e-5; //Dna junc to SL [cm^2/sec]
        const double D_na_SLcyto = 1.79e-5; //Dna SL to cyto [cm^2/sec]
        const double J_ca_juncsl = (D_ca_JuncSL * SAjunc / dist_JuncSL * 1e-10); //[L/msec] [m^2/sec] MOUSE
        const double J_ca_slmyo = (DcaSLcyto * SAsl / distSLcyto * 1e-10);
        const double J_na_juncsl = (D_na_JuncSL * SAjunc / dist_JuncSL * 1e-10);
        const double J_na_slmyo = (D_na_SLcyto * SAsl / distSLcyto * 1e-10);

        //Temperature correction
        const double Qpow = ((sim_config.temp - 295.0) / 10.0);
        const double Q10to_timeconstant = 1.8; //1.4+-0.2<<Brouillette,2004/, 2.66<<Noujaim,2007/, 2.1<<Courtemanche,1998/, 2.1<<(Benndorf and Nilius, 1987;Schwarz, 1986/ !!!!
        const double Q10to_conductance = 1.1;
        const double Q10Kur_timeconstant = 1.9; //1.5+-0.3<<Brouillette,2004/, 2.2<<Courtemanche,1998/ !!!! by shanzhuo
        const double Q10Kur_conductance = 1.37;
        const double Q10Kss = 1.9; //1.9+-0.4<<Brouillette,2004/ !!!!
        //const double Q10CaL = 1; //not in use !!!!
        const double Q10_tau_h_Na = 1.0 / 2.5;

        // Fixed ion concentrations
        const double _Nao_ = 140.0; //Extracellular Na  [mM]
        const double _Ko_  = 5.4; //Extracellular K   [mM]
        const double _Cao_ = 1.0; //Extracellular Ca  [mM] MOUSE // 1.8 mM in RABBIT
        const double _Clo_ = 145.0; //Extracellular Cl  [mM]
        const double _Cli_ = 15.0; //Intracellular Cl  [mM]
        const double _Mgi_ = 1.0; //Intracellular Mg  [mM]

        //Conductance:

        double g_Kto = (1 * 0.2320);    //[nS/pF] Ito=11.2  (+25mV, 22°C) <<Lomax,2003>> !!!!
        double g_Kur = (1 * 0.09795);   //IKur=4.68 (+20mV, 22°C) <<TB, 2004>> !!!!
        double g_Kss = (1 * 0.049);     //IKss=3.67 (+30mV, 22°C) <<Xu, 1999>> !!!!
        double g_Kr = (1 * 0.06328);    //IKr=2.5   (+15mV, 37°C) <<Nakamura, 2010>> !!!!
        double g_KACh = (6.534);        //
        double g_K1 = (1 * 0.5291);     // (22°C), <<Lomax,2003>>, right shifted 3mV at 37°C !!!!

        double ICa_scale_ratio = (0.5 * 0.5 * 1);// !!!! (weijian & shanzhuo)
        double _IbarNCX_ = (1.0 * 1);
        double _IbarNaK_ = (0.75 * 0.8 * 5.0);// !!!! (weijian & shanzhuo)
        double _IbarpCa_ = (0.8 * 0.0673);// !!!! (weijian)
        double _GNa_ = (10);
        double _GCaB_ = (1 * 7.539e-4);
        double _GNaB_ = (1 * 4.0 * 0.297e-3);
        double _GKB_ = 0.01;

        //Ca2+ flux:

        //Ca transport parameters
        double _KmCai_ = (3.59e-3);   // [mM]
        double _KmCao_ = 1.3; // [mM]
        double _KmNai_ = (12.29);   // [mM]
        double _KmNao_ = 87.5; // [mM]

        double _ksat_ = (0.27);   // [none]

        double _nu_ = 0.35; // [none]
        double _Kdact_ = 0.128e-3; // [mM] changed from rabbit 1.0/2.0*0.256e-3
        double _KmPCa_ = 0.5e-3; // [mM]

        //SERCA
        double _hillSRCaP_ = 1.787;     // [mM]
        double _Kmr_ = 2.1;      // [mM] changed from rabbit (1.7) // from Yang-Saucerman
        double _Kmf_ = (0.337e-3 + 0 * 0.3e-3);   // [mM] changed from rabbit (0.246e-3) // from Yang-Saucerman !!!!
        double _Vmax_SERCA_ = 0.5 * (1.15 * 1.15 * 0.000286);// [mM/msec] (mmol/L cytosol/msec) // Morotti // * 1/2 Shanzhuo

        //RyR
        double _koCa_ = 10.0; // [mM^-2 1/ms]
        double _kiCa_ = 0.5; // [1/mM/ms]
        double _kom_ = 0.06; // [1/ms]
        double _kim_ = 0.005; // [1/ms]
        double _ec50SR_ = (0.45 + 0 * 0.50); // [mM] changed from rabbit (0.45) !!!! weijian
        double _ks_ = 25.0; // [1/ms]

        //Signaling Pathways Parameters:

        const double _plb_val_ = 106.0; //RABBIT: 38.0
        //Parameters for CaMKII module
        const double _LCCtotJunc_ = 31.4*0.9; //[um] - Total Juncic [LCC] - (umol/l junc)
        const double _LCCtotSL_ = 0.0846; //[um] - Total Subsarcolemmal [LCC] (umol/l sl)
        const double _RyRtot_ = 382.6; //[um] - Total RyR (in Junc)
        const double _PP1_junc_ = 95.7; //[um] - Total junctional [PP1]
        const double _PP1_SL_ = 0.57; //[um] - Total Subsarcolemmal [PP1]
        const double _PP2A_junc_ = 95.76; //[um] - Total junctional PP2A
        const double _OA_ = 0; //[um] - PP1/PP2A inhibitor Okadaic Acid
        const double _PLBtot_ = _plb_val_; //[um] - Total [PLB] in cytosolic units
 
        //Parameters for BAR module
        const double _LCCtotBA_ = 0.025; //[um] - [umol/L cytosol]
        const double _RyRtotBA_ = 0.135; //[um] - [umol/L cytosol]
        const double _PLBtotBA_ = _plb_val_; //[um] - [umol/L cytosol]
        const double _TnItotBA_ = 70.0; //[um] - [umol/L cytosol]
        const double _IKstotBA_ = 0.025; //[um] - [umol/L cytosol]
        const double _ICFTRtotBA_ = 0.025; //[um] - [umol/L cytosol]
        const double _PP1_PLBtot_ = 0.89; //[um] - [umol/L cytosol]
        const double _IKurtotBA_ = 0.025; //[um] - [umol/L cytosol] MOUSE
        const double _PLMtotBA_ = 48.0; //[um] - [umol/L cytosol] MOUSE

        //CaMKII phos Parameters (all units in um)
        const double _BtotJunc_ = 1.54 / 8.293e-4;
        const double _kSLmyo_ = 8.587e-15; //[L/msec]
        const double _CaNtotJunc_ = 3e-3 / 8.293e-4;
        const double _CaNtotSL_ = 3e-3;
        const double _CaNtotCyt_ = 3e-3;
        const double _PP1totJunc_ = 96.5;
        const double _PP1totSL_ = 0.57;
        const double _PP1totCyt_  = 0.57;
        const double _BtotSL_ = 24.2;
        const double _BtotCyt_ = 24.2;
        const double _n_OE_ = 6.0;
        double _CaMKIItotJunc_; //Initialised within CaM_ODEs, dependent on "_CaMKII_level_" config
        double _CaMKIItotSL_; //As above
        double _CaMKIItotCyt_; //As above

        //CaM buffering (B) parameters
        const double _k0Boff_ = 0.0014; //CaM dissociation from CaMB (Buffer) [s^-1]
        const double _k0Bon_ = _k0Boff_/0.2; //CaM association to CaMB (Buffer) [um^-1 s^-1] kon = koff/Kd
        const double _k2Boff_ = _k0Boff_/100.0; //Ca2CaM dissocation from Ca2CaMB [s^-1]
        const double _k2Bon_ = _k0Bon_; //Ca2CaM association to Ca2CaMB [um^-1 s^-1]
        const double _k4Boff_ = _k2Boff_; //Ca4CaM dissocation from Ca4CaMB [s^-1]
        const double _k4Bon_ = _k0Bon_; //Ca4CaM association to Ca4CaMB [um^-1 s^-1]

        //Buffer Parameters:

        //TnC low affinity
        const double _Bmax_TnClow_ = 70e-3; //[mM]
        const double _koff_tncl_ = 19.6e-3; //[1/ms]
        const double _kon_tncl_ = 32.7; //[1/mM/ms]
        //TnC high affinity
        const double _Bmax_TnChigh_ = 140e-3; //[mM]
        const double _koff_tnchca_ = 0.032e-3; //[1/ms]
        const double _kon_tnchca_ = 2.37; //[1/mM/ms]
        const double _koff_tnchmg_ = 3.33e-3; //[1/ms]
        const double _kon_tnchmg_ = 3e-3; //[1/mM/ms]
        //Myosin buffering
        const double _Bmax_myosin_ = 140e-3; //[mM]
        const double _koff_myoca_ = 0.46e-3; //[1/ms]
        const double _kon_myoca_ = 13.8; //[1/mM/ms]
        const double _koff_myomg_ = 0.057e-3; //[1/ms]
        const double _kon_myomg_ = 0.0157; //[1/mM/ms]
        //SRB buffering
        const double _Bmax_SR_ = 17.1e-3; //[mM]
        const double _koff_sr_ = 60e-3; //[1/ms]
        const double _kon_sr_ = 100.0; //[1/mM/ms]
        //CSQN
        const double _Bmax_Csqn_ = 2.7; //[mM]
        const double _koff_csqn_ = 65.0; //[1/ms]
        const double _kon_csqn_ = 100.0; //[1/mM/ms]
        //Na buffering
        const double _Bmax_Naj_ = 7.561; //[mM]
        const double _Bmax_Nasl_ = 1.65; //[mM]
        const double _koff_na_ = 1e-3; //[1/ms]
        const double _kon_na_ = 0.1e-3; //[1/mM/ms]
        //SL buffering Low
        const double _Bmax_SLlowsl_ = 1.214845; //(37.38e-3*Vmyo/Vsl)       atria:1.74   // [mM]
        const double _Bmax_SLlowj_ = 0.557143; //(4.62e-3*Vmyo/Vjunc*0.1)  atria:1.83    // [mM]
        const double _koff_sll_ = 1300e-3; //[1/ms]
        const double _kon_sll_ = 100.0; //[1/mM/ms]
        //SL buffering High
        const double _Bmax_SLhighsl_ = 0.433875; //(13.35e-3*Vmyo/Vsl)        atria:0.62   // [mM]
        const double _Bmax_SLhighj_ = 0.198980; //(1.65e-3*Vmyo/Vjunc*0.1)  atria:0.65 // [mM]
        const double _koff_slh_ = 30e-3; //[1/ms]
        const double _kon_slh_ = 100.0; //[1/mM/ms]
        //Caculations
        const double KQ10_to = pow(Q10to_timeconstant,Qpow);
        const double KQ10_Kur = pow(Q10Kur_timeconstant,Qpow);
        const double KQ10_Kss = pow(Q10Kss,Qpow);

        //Simulation parameters
        double dVdt_max{-5000}; //Holds max dV/dt

        double V_old{}; //Potential from prev timestep
        double dVdt_old{}; //Voltage rate for the previous timestep, for use in measurements
        double cycle_start; 

        double MDP{}; //Resting membrane potential
        double min_V{}; //Significant potentials over coourse of AP
        double max_V{};
        double APD_start{}; //APD start times for a given AP
        double APD90_V{}; //Calculated voltages for each of the APDs
        double APD50_V{};
        double APD30_V{};
        double APD30{};
        double APD50{};
        double APD90{};
        double norm_APD50{};

        double I_j; // [pA/pF] Junction current
        
        //Declare output files that use fprintf for formatting
        FILE *output_currents_file;
        FILE *output_measurements_file;

    public:
        atrial_cell() {}
        ~atrial_cell() {}
        const double &get_potential() {return states.ec.V;}
        void set_I_j(const double &junc_current) 
        {I_j = junc_current/(Cmem*1e12);} //[->A/F]
        void setup_simulation();
        void output_config(std::ofstream &);
        void apply_config();
        void create_data_file();
        void create_measurement_file();
        void set_initial_states();
        void store_variables();
        void update_state_variables();
        void ODEs(const double &);
        void IP3_ODEs(const double &);
        void cam_ODEs(const double &);
        void cam_compart(const double *, const cams *, cams *);
        void camkii_ODEs(const double &);
        void bars_ODEs(const double &);
        void ecc_ODEs(const double &);
        double INa_channel(const double &);
        double INaK_channel(const double &);
        double Ito_channel(const double &);
        double IKur_channel();
        double IKss_channel();
        double IKr_channel();
        double IKACh_channel();
        double IK1_channel();
        double IKCa_channel();
        double background_currents(const double &, const int);
        double INCX_channel();
        double LCC_channel();
        void LCC_rates(const double *, double *);
        void Na_concentrations();
        void classic_Ca_handling(const double &);
        void IP3R_rates(const double *, double *);
        double RyR_rel(const double &, const double *, double *);
        double RyR_leak(const double *);
        double SERCA(const double &, const double *);
        void Euler_method();
        void output_data(const double &);
        void output_final_states();
        void measurements(const double &, int &);
        void output_measurements(const int &);
};

#endif