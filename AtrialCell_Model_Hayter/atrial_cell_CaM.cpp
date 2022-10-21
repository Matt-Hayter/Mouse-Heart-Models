#include "atrial_cell.hpp"

//ODEs for cams
void atrial_cell::cam_ODEs()
{
    //Perform shallow copy of relevant state variables and rates of change, for conciseness
    const eccs *ec = &(states.ec);
    const cams *cj = &(states.cj);
    const cams *cs = &(states.cs);
    const cams *cc = &(states.cc);
    cams *cjR = &(states_dot.cj);
    cams *csR = &(states_dot.cs);
    cams *ccR = &(states_dot.cc);

    double CaMKIItotJunc = (config.ACUTE == 1 && time > config.ACUTE_TIME) ? _n_OE_ * _CaMKIItotJunc_ : _CaMKIItotJunc_;
    double pCaMCyt[] = {_BtotCyt_, _CaMKIItotCyt_, _CaNtotCyt_, _PP1totCyt_, states.ec.Ca_i * 1e3, 0};
    double pCaMSL[] = {_BtotSL_, _CaMKIItotSL_, _CaNtotSL_, _PP1totSL_, states.ec.Ca_sl * 1e3, 1};
    double pCaMJunc[] = {0, _CaMKIItotJunc_, _CaNtotJunc_, _PP1totJunc_, states.ec.Ca_junc * 1e3, 2};

    //Process ODEs for cc, cs and cj cam states individually
    cam_compart(pCaMCyt, cc, ccR);
    cam_compart(pCaMSL, cs, csR);
    cam_compart(pCaMJunc, cj, cjR);

    //// CaM diffusion between modules
    double CaMtotJunc_o = (cj->CaM + cj->Ca2CaM + cj->Ca4CaM + cj->CaMB + cj->Ca2CaMB + cj->Ca4CaMB)
        + (cj->Pb2 + cj->Pb + cj->Pt + cj->Pt2) * CaMKIItotJunc + (cj->CaMCa4CaN + cj->Ca4CaMCa4CaN);

    double Bjunc = _BtotJunc_ - CaMtotJunc_o; // [uM junc]
    double J_cam_juncSL = 1e-3 * (_k0Boff_ * cj->CaM - _k0Bon_ * Bjunc * cs->CaM); // [uM/msec junc]
    double J_ca2cam_juncSL = 1e-3 * (_k2Boff_ * cj->Ca2CaM - _k2Bon_ * Bjunc * cs->Ca2CaM); // [uM/msec junc]
    double J_ca4cam_juncSL = 1e-3 * (_k4Boff_ * cj->Ca4CaM - _k4Bon_ * Bjunc * cs->Ca4CaM); // [uM/msec junc]
    double J_cam_SLmyo = _kSLmyo_ * (cs->CaM - cc->CaM); // [umol/msec]
    double J_ca2cam_SLmyo = _kSLmyo_ * (cs->Ca2CaM - cc->Ca2CaM); // [umol/msec]
    double J_ca4cam_SLmyo = _kSLmyo_ * (cs->Ca4CaM - cc->Ca4CaM); // [umol/msec]

    cjR->CaM -= J_cam_juncSL;
    cjR->Ca2CaM -= J_ca2cam_juncSL;
    cjR->Ca4CaM -= J_ca4cam_juncSL;
    csR->CaM += J_cam_juncSL * Vjunc / Vsl - J_cam_SLmyo / Vsl;
    csR->Ca2CaM += J_ca2cam_juncSL * Vjunc / Vsl - J_ca2cam_SLmyo / Vsl;
    csR->Ca4CaM += J_ca4cam_juncSL * Vjunc / Vsl - J_ca4cam_SLmyo / Vsl;
    ccR->CaM += J_cam_SLmyo / Vmyo;
    ccR->Ca2CaM += J_ca2cam_SLmyo / Vmyo;
    ccR->Ca4CaM += J_ca4cam_SLmyo / Vmyo;
}

//ODEs are process for all cams passed to this function
void atrial_cell::cam_compart(const double *p, const cams *cam, cams *camR)
{
    const eccs *ec = &(states.ec);
    eccs *ecR = &(states_dot.ec);

    const double Ki = ec->Ki;

    const double Btot = p[0];
    const double CaMKIItot = p[1];
    const double CaNtot = p[2];
    const double PP1tot = p[3];
    const double Ca = p[4]; //[uM]
    const int type = p[5];

    // ---------- Ca/CaM Module ----------
    double Kd02, Kd24;
    if (_Mgi_ > 1.0) {
        Kd02 = 0.0025 * (1 + Ki / 0.94 - 1 / 0.012 + (_Mgi_ - 1) / 0.060) * (1 + Ki / 8.1 + 1 / 0.022 + (_Mgi_ - 1) / 0.068);   // [uM^2]
        Kd24 = 0.128 * (1 + Ki / 0.64 + 1 / 0.0014 + (_Mgi_ - 1) / 0.005) * (1 + Ki / 13.0 - 1 / 0.153 + (_Mgi_ - 1) / 0.150);  // [uM^2]
    } else {
        Kd02 = 0.0025 * (1 + Ki / 0.94 - _Mgi_ / 0.012) * (1 + Ki / 8.1 + _Mgi_ / 0.022);  // [uM^2]
        Kd24 = 0.128 * (1 + Ki / 0.64 + _Mgi_ / 0.0014) * (1 + Ki / 13.0 - _Mgi_ / 0.153); // [uM^2]
    }
    double k20 = 10;               // 2 Ca dissociation from Ca2CaM     [s^-1]
    double k02 = k20 / Kd02;   // 2 Ca association with CaM         [uM^-2 s^-1]
    double k42 = 500;              // 2 Ca dissociation from Ca4CaM     [s^-1]
    double k24 = k42 / Kd24;   // 2 Ca association with Ca2CaM      [uM^-2 s^-1]

    //Using thermodynamic constraints
    double k20B = k20 / 100; // 2 Ca dissociation from Ca2CaMB (Buffer)     [s^-1] thermo constraint on loop 1
    double k02B = k02; // 2 Ca association with CaMB                  [uM^-2 s^-1]
    double k42B = k42; // 2 Ca dissociation from Ca4CaMB              [s^-1] thermo constraint on loop 2
    double k24B = k24; // 2 Ca association with Ca2CaMB               [uM^-2 s^-1]

    // CaM Reaction fluxes
    double rcn02 = k02 * pow(Ca, 2) * cam->CaM - k20 * cam->Ca2CaM;
    double rcn24 = k24 * pow(Ca, 2) * cam->Ca2CaM - k42 * cam->Ca4CaM;

    // CaM buffer fluxes
    double B = Btot - cam->CaMB - cam->Ca2CaMB - cam->Ca4CaMB;
    double rcn02B = k02B * pow(Ca, 2) * cam->CaMB - k20B * cam->Ca2CaMB;
    double rcn24B = k24B * pow(Ca, 2) * cam->Ca2CaMB - k42B * cam->Ca4CaMB;
    double rcn0B = _k0Bon_ * cam->CaM * B - _k0Boff_ * cam->CaMB;
    double rcn2B = _k2Bon_ * cam->Ca2CaM * B - _k2Boff_ * cam->Ca2CaMB;
    double rcn4B = _k4Bon_ * cam->Ca4CaM * B - _k4Boff_ * cam->Ca4CaMB;

    // ---------- CaMKII Module ----------
    double kbi = 2.2;               // Ca4CaM dissocation Pb -> Pi     [s^-1]
    double kib = kbi / 33.5e-3;  // Ca4CaM assocation  Pi -> Pb2    [uM^-1 s^-1]
    double kib2 = kib;           // Ca2CaM dissocation Pb2 -> Pi    [uM^-1 s^-1]
    double kb2i = kib2 * 5;           // Ca2CaM assocation  Pi -> Pb2    [s^-1]
    double kb24 = k24;           // 2 Ca association   Pb2 -> Pb    [uM^-2 s^-1]
    double kb42 = k42 * 33.5e-3 / 5.0; // 2 Ca dissociation  Pb -> Pb2    [s^-1]
    double kpp1 = 1.72;       // PP1-dep (Thr287) dephosphorylation rates [s^-1]
    double Kmpp1 = 11.5;        // PP1-dep (Thr287) dephosphorylation rates [uM]
    double kta = kbi / 1000;          // Ca4CaM dissociation Pt -> Pa    [s^-1]
    double kat = kib;            // Ca4CaM association  Pa -> Pt    [uM^-1 s^-1]
    double kt42 = k42 * 33.5e-6 / 5.0; // 2 Ca dissociation   Pt -> Pt2   [s^-1]
    double kt24 = k24;           // 2 Ca association    Pt2 -> Pt   [uM^-2 s^-1]
    double kat2 = kib;           // Ca2CaM association  Pa -> Pt2   [uM^-1 s^-1]
    double kt2a = kib * 5;            // Ca2CaM dissociation Pt2 -> Pa   [s^-1]

    // CaMKII reaction fluxes
    double Pi = 1 - cam->Pb2 - cam->Pb - cam->Pt - cam->Pt2 - cam->Pa;
    double rcnCKib2 = kib2 * cam->Ca2CaM * Pi - kb2i * cam->Pb2;
    double rcnCKb2b = kb24 * pow(Ca, 2) * cam->Pb2 - kb42 * cam->Pb;
    double rcnCKib = kib * cam->Ca4CaM * Pi - kbi * cam->Pb;
    double T = cam->Pb + cam->Pt + cam->Pt2 + cam->Pa;
    double kbt = 0.055 * T + .0074 * pow(T, 2) + 0.015 * pow(T, 3);
    double rcnCKbt = kbt * cam->Pb - kpp1 * PP1tot * cam->Pt / (Kmpp1 + CaMKIItot * cam->Pt);
    double rcnCKtt2 = kt42 * cam->Pt - kt24 * pow(Ca, 2) * cam->Pt2;
    double rcnCKta = kta * cam->Pt - kat * cam->Ca4CaM * cam->Pa;
    double rcnCKt2a = kt2a * cam->Pt2 - kat2 * cam->Ca2CaM * cam->Pa;
    double rcnCKt2b2 = kpp1 * PP1tot * cam->Pt2 / (Kmpp1 + CaMKIItot * cam->Pt2);
    double rcnCKai = kpp1 * PP1tot * cam->Pa / (Kmpp1 + CaMKIItot * cam->Pa);

    // ---------- CaN Module ----------
    double kcanCaoff = 1;    // 2 Ca dissociation from Ca4CaN -> Ca2CaN   [s^-1]
    double kcanCaon = kcanCaoff / 0.5; // 2 Ca association  from Ca2CaN -> Ca4CaN   [uM^-1 s^-1]
    double kcanCaM4on = 46;            // [uM^-1 s^-1]
    double kcanCaM4off = 1.3e-3;       // [s^-1]
    double kcanCaM2on = kcanCaM4on;
    double kcanCaM2off = 2508 * kcanCaM4off;
    double kcanCaM0on = kcanCaM4on;
    double kcanCaM0off = 165 * kcanCaM2off;
    double k02can = k02;
    double k20can = k20 / 165.0;
    double k24can = k24;
    double k42can = k20 / 2508.0;

    // CaN reaction fluxes
    double Ca2CaN = CaNtot - cam->Ca4CaN - cam->CaMCa4CaN - cam->Ca2CaMCa4CaN - cam->Ca4CaMCa4CaN;
    double rcnCa4CaN = kcanCaon * pow(Ca, 2) * Ca2CaN - kcanCaoff * cam->Ca4CaN;
    double rcn02CaN = k02can * pow(Ca, 2) * cam->CaMCa4CaN - k20can * cam->Ca2CaMCa4CaN;
    double rcn24CaN = k24can * pow(Ca, 2) * cam->Ca2CaMCa4CaN - k42can * cam->Ca4CaMCa4CaN;
    double rcn0CaN = kcanCaM0on * cam->CaM * cam->Ca4CaN - kcanCaM0off * cam->CaMCa4CaN;
    double rcn2CaN = kcanCaM2on * cam->Ca2CaM * cam->Ca4CaN - kcanCaM2off * cam->Ca2CaMCa4CaN;
    double rcn4CaN = kcanCaM4on * cam->Ca4CaM * cam->Ca4CaN - kcanCaM4off * cam->Ca4CaMCa4CaN;

    //// ---------- Equations ----------

    // CaM equations
    double dCaM = 1e-3 * (-rcn02 - rcn0B - rcn0CaN);
    double dCa2CaM = 1e-3 * (rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItot * (-rcnCKib2 + rcnCKt2a));
    double dCa4CaM = 1e-3 * (rcn24 - rcn4B - rcn4CaN + CaMKIItot * (-rcnCKib + rcnCKta));
    double dCaMB = 1e-3 * (rcn0B - rcn02B);
    double dCa2CaMB = 1e-3 * (rcn02B + rcn2B - rcn24B);
    double dCa4CaMB = 1e-3 * (rcn24B + rcn4B);

    // CaMKII equations
    double dPb2 = 1e-3 * (rcnCKib2 - rcnCKb2b + rcnCKt2b2);
    double dPb = 1e-3 * (rcnCKib + rcnCKb2b - rcnCKbt);
    double dPt = 1e-3 * (rcnCKbt - rcnCKta - rcnCKtt2);
    double dPt2 = 1e-3 * (rcnCKtt2 - rcnCKt2a - rcnCKt2b2);
    double dPa = 1e-3 * (rcnCKta + rcnCKt2a - rcnCKai);

    // CaN equations
    double dCa4CaN = 1e-3 * (rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);
    double dCaMCa4CaN = 1e-3 * (rcn0CaN - rcn02CaN);
    double dCa2CaMCa4CaN = 1e-3 * (rcn2CaN + rcn02CaN - rcn24CaN);
    double dCa4CaMCa4CaN = 1e-3 * (rcn4CaN + rcn24CaN);

    camR->CaM = dCaM;
    camR->Ca2CaM = dCa2CaM;
    camR->Ca4CaM = dCa4CaM;
    camR->CaMB = dCaMB;
    camR->Ca2CaMB = dCa2CaMB;
    camR->Ca4CaMB = dCa4CaMB;
    camR->Pb2 = dPb2;
    camR->Pb = dPb;
    camR->Pt = dPt;
    camR->Pt2 = dPt2;
    camR->Pa = dPa;
    camR->Ca4CaN = dCa4CaN;
    camR->CaMCa4CaN = dCaMCa4CaN;
    camR->Ca2CaMCa4CaN = dCa2CaMCa4CaN;
    camR->Ca4CaMCa4CaN = dCa4CaMCa4CaN;

    // write to global variables for adjusting Ca buffering in EC coupling model
    double JCa = 1e-3 * 1e-3 * 2 * (CaMKIItot * (rcnCKtt2 - rcnCKb2b) - (rcn02 + rcn24 + rcn02B + rcn24B + rcnCa4CaN + rcn02CaN + rcn24CaN)); // [mM/msec]
    if (type == 0) {
        ecR->Ca_i += JCa; //JCaCyt
    } else if (type == 1) {
        ecR->Ca_sl += JCa; //JCaSL
    } else {
        ecR->Ca_junc += JCa; //JCaJunc
    }
}