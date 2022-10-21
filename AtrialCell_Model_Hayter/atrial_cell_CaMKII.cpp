#include "atrial_cell.hpp"

//ODEs for camkiis
void atrial_cell::camkii_ODEs()
{
    //Perform shallow copy of relevant state variables and rates of change, for conciseness
    const cams *cj = &(states.cj);
    const cams *cs = &(states.cs);
    const cams *cc = &(states.cc);
    const bars *ba = &(states.ba);
    const ckiis *ck = &(states.ck);
    ckiis *ckR = &(states_dot.ck);
    
    double CaMKIItotJunc = (config.ACUTE == 1 && time > config.ACUTE_TIME) ? _n_OE_ * _CaMKIItotJunc_ : _CaMKIItotJunc_;
    double CaMKIItotSL = (config.ACUTE == 1 && time > config.ACUTE_TIME) ? _n_OE_ * _CaMKIItotSL_ : _CaMKIItotSL_;
    double CaMKIItotCyt = (config.ACUTE == 1 && time > config.ACUTE_TIME) ? _n_OE_ * _CaMKIItotCyt_ : _CaMKIItotCyt_;

    double CaMKIIact_Junc = CaMKIItotJunc * (cj->Pb + cj->Pt + cj->Pt2 + cj->Pa); // Multiply total by fraction
    double CaMKIIact_SL = CaMKIItotSL * (cs->Pb + cs->Pt + cs->Pt2 + cs->Pa);
    double CaMKIIact_Cyt = CaMKIItotCyt * (cc->Pb + cc->Pt + cc->Pt2 + cc->Pa);
    double PP1_PLB_avail = 1 - ba->I1p_PP1 / _PP1_PLBtot_ + 0.081698; // NEW ODEs

    dependents.camkii_junc_act = CaMKIIact_Junc;
    dependents.camkii_cyt_act = CaMKIIact_Cyt;
    dependents.camkii_sl_act = CaMKIIact_SL;

    // double LCC_PKAp		= ck->LCC_PKAp		;
    double LCC_CKjuncp = ck->LCC_CKjuncp; // LCCp-CaMKIIjunc , Juncic [LCCp] by junctional CaMKII
    // double RyR2809p		= ck->RyR2809p		;	// RyR-Ser2809p    , [RyR-Ser2809p] by PKA (currently unused anywhere else)
    double RyR2815p = ck->RyR2815p;    // RyR-Ser2815p    , [RyR-Ser2815p] by CaMKII
    double PLBT17p = ck->PLBT17p;    // PLB-Thr17p      , [PLB-Thr17p] by CaMKII
    double LCC_CKslp = ck->LCC_CKslp;// LCCp-CaMKIIsl   , Subsarcolemmal [LCCp] by subsarcolemmal CaMKII

    // L-Type Ca Channel (LTCC) parameters
    double k_ckLCC = 0.4;         // [s^-1]	CaMKII phosphorylation rate at LCC
    double k_pp1LCC = 0.1103;    // [s^-1] 	PP1 dephosphorylation rate at LCC
    //double k_pkaLCC   = 13.5;                // [s^-1]
    //double k_pp2aLCC  = 10.1;               // [s^-1]
    double KmCK_LCC = 12.0; // [uM] 		Michaelis Constant for CaMKII phosphorylation at LCC
    //double KmPKA_LCC  = 21.0;               // [uM]
    //double KmPP2A_LCC = 47.0;              // [uM]
    double KmPP1_LCC = 9.0; // [uM] 		Michaelis Constant for PP1 dephosphorylation at LCC

    // Ryanodine Receptor (RyR) parameters
    double k_ckRyR = 0.4;        // [s^-1] 	CaMKII phosphorylatino rate at RyR
    //double k_pkaRyR  = 1.35;                // [s^-1]
    double k_pp1RyR = 1.07;      // [s^-1] 	PP1 dephosphorylation rate at RyR
    double k_pp2aRyR = 0.481;    // [s^-1] 	PP2A dephosphorylation rate at RyR

    // Basal RyR phosphorylation (numbers based on param estimation)
    //double kb_2809    = 0.51;                 // [uM/s] 	PKA site
    double kb_2815 = 0.35; // [uM/s] 	CaMKII site, Basal phosphorylation of Ser2815
    double KmCK_RyR = 12.0; // [uM] 		Michaelis Constant for CaMKII phosphorylation at RyR
    //double KmPKA_RyR  = 21.0;               // [uM]
    double KmPP1_RyR = 9.0; // [uM] 		Michaelis Constant for PP1 phosphorylation at RyR
    double KmPP2A_RyR = 47.0; // [uM] 		Michaelis Constant for PP2A phosphorylation at RyR

    //Phospholamban (PLB) parameters
    double k_ckPLB = 8e-3;        // [s^-1]	CaMKII phosphorylation rate at PLB
    double k_pp1PLB = 0.0428;     // [s^-1]	PP1 dephosphorylation rate at PLB
    double KmCK_PLB = 12.0;    // [uM]		Michaelis Constant for CaMKII phosphorylation at PLB
    double KmPP1_PLB = 9.0;    // [uM]		Michaelis Constant for PP1 phosphorylation at PLB

    //Okadaic Acid inhibition params (based on Huke/Bers [2008]), treat _OA_ as non-competitive inhibitor of PP1 and PP2A
    double Ki_OA_PP1 = 0.78; // [uM] - Values from fit, Inhibitory constant of _OA_ on PP1
    double Ki_OA_PP2A = 0.037; // [uM] - Values from fit, Inhibitory constant of _OA_ on PP2A

    // Default PKA level
    //double PKAc = 95.6 * .54;

    //_OA_ inhibition term (non-competitive) for PP1 and PP2A
    double OA_PP1 = 1 / (1 + pow(_OA_ / Ki_OA_PP1, 3));
    double OA_PP2A = 1 / (1 + pow(_OA_ / Ki_OA_PP2A, 3));

    //ODE equations:

    //LTCC states (note: PP2A is acting on PKA site and PP1 on CKII site)
    //CaMKII phosphorylation of Juncic LCCs
    double LCC_CKjuncn = _LCCtotJunc_ - LCC_CKjuncp;
    double LCCJunc_PHOS = (k_ckLCC * CaMKIIact_Junc * LCC_CKjuncn) / (KmCK_LCC + LCC_CKjuncn);
    double LCCJunc_DEPHOS = (k_pp1LCC * _PP1_junc_ * LCC_CKjuncp) / (KmPP1_LCC + LCC_CKjuncp) * OA_PP1;
    double dLCC_CKjuncp = LCCJunc_PHOS - LCCJunc_DEPHOS;

    //CaMKII phosphorylation of Sub-sarcolemmal LCCs
    double LCC_CKsln = _LCCtotSL_ - LCC_CKslp;
    double LCCSL_PHOS = (k_ckLCC * CaMKIIact_SL * LCC_CKsln) / (KmCK_LCC + LCC_CKsln);
    double LCCSL_DEPHOS = (k_pp1LCC * _PP1_SL_ * LCC_CKslp) / (KmPP1_LCC + LCC_CKslp) * OA_PP1;
    double dLCC_CKslp = LCCSL_PHOS - LCCSL_DEPHOS;

    //RyR states
    double RyR2815n = _RyRtot_ - RyR2815p;
    double RyR_BASAL = kb_2815 * RyR2815n;
    double RyR_PHOS = (k_ckRyR * CaMKIIact_Junc * RyR2815n) / (KmCK_RyR + RyR2815n);
    double RyR_PP1_DEPHOS = (k_pp1RyR * _PP1_junc_ * RyR2815p) / (KmPP1_RyR + RyR2815p) * OA_PP1;
    double RyR_PP2A_DEPHOS = (k_pp2aRyR * _PP2A_junc_ * RyR2815p) / (KmPP2A_RyR + RyR2815p) * OA_PP2A;
    double dRyR2815p = RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS;

    //// PLB states
    double PP1_PLB = _PP1_junc_ * PP1_PLB_avail; // Inhibitor-1 regulation of _PP1_junc_ included here
    double PLBT17n = _PLBtot_ - PLBT17p;
    double PLB_PHOS = (k_ckPLB * PLBT17n * CaMKIIact_Junc) / (KmCK_PLB + PLBT17n);
    double PLB_DEPHOS = (k_pp1PLB * PP1_PLB * PLBT17p) / (KmPP1_PLB + PLBT17p) * OA_PP1;
    double dPLBT17p = PLB_PHOS - PLB_DEPHOS;

    //// Collect ODEs and convert to uM/ms
    // ckR->LCC_PKAp	 = 1e-3 * dLCC_PKAp;
    ckR->LCC_CKjuncp = 1e-3 * dLCC_CKjuncp;
    // ckR->RyR2809p	 = 1e-3 * dRyR2809p;
    ckR->RyR2815p = 1e-3 * dRyR2815p;
    ckR->PLBT17p = 1e-3 * dPLBT17p;
    ckR->LCC_CKslp = 1e-3 * dLCC_CKslp;

    //// CaMKII phos
    dependents.ph.LCC_CKp = ck->LCC_CKjuncp / _LCCtotJunc_; // fractional CaMKII-dependent LCC junc phos
    dependents.ph.RyR_CKp = ck->RyR2815p / _RyRtot_; // fractional CaMKII-dependent RyR phos
    dependents.ph.PLB_CKp = ck->PLBT17p / _PLBtot_; // fractional CaMKII-dependent PLB phos
}