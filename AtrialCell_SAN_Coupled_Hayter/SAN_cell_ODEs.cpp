#include "SAN_cell.hpp"
#include <cmath>

//SAN ODEs
void SAN_cell::ODEs(const double &time)
{
    //For conciseness
    const double T = sim_config.temp;
    //Copy over state variables for conciseness
    double v = states.v;
    double dst = states.dst;
    double fst = states.fst;
    double dt  = states.dt ;
    double ft = states.ft;
    double ikr_act;
    double ikr_act_f = states.ikr_act_f;
    double ikr_act_s = states.ikr_act_s;
    double ikr_inact = states.ikr_inact;
    double iks_act = states.iks_act;
    double fl12 = states.fl12;
    double dl12 = states.dl12;
    double fl13 = states.fl13;
    double dl13 = states.dl13;
    double r = states.r;
    double m_ttxr = states.m_ttxr;
    double h_ttxr = states.h_ttxr;
    double j_ttxr = states.j_ttxr;
    double m_ttxs = states.m_ttxs;
    double h_ttxs = states.h_ttxs;
    double j_ttxs = states.j_ttxs;
    double y_1 = states.y_1;
    double carel = states.carel;
    double caup = states.caup;
    double casub = states.casub;
    double Ftc = states.Ftc;
    double Ftmc = states.Ftmc;
    double Ftmm = states.Ftmm;
    double Fcms = states.Fcms;
    double Fcmi = states.Fcmi;
    double Fcq = states.Fcq;
    double cai = states.cai;
    double q = states.q;
    double fca = states.fca;
    double nai = states.nai;
    double ki = states.ki;
    double resting = states.resting;
    double open = states.open;
    double inactivated = states.inactivated;
    double resting_inactivated = states.resting_inactivated;
    double y_2 = states.y_2;
    double y_4 = states.y_4;

    /*Reversal Potentials****************************************************/

        double ena,eca,ek,eks;

        ena = (R*T/F)*log(nao/nai);
        ek  = (R*T/F)*log(ko/ki);
        eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
        eca = (R*T/(2*F))*log(cao/casub);

    /*Ist********************************************************************/

        double qa, qi, tauqa,tauqi,alphaqa,betaqa,alphaqi,betaqi;
        qa = 1.0/(1.0 + exp(-(v+67.0)/5.0));
        alphaqa = 1.0/(0.15*exp(-(v)/11.0)+0.2*exp(-(v)/700.0));
        betaqa  =  1.0/(16.0*exp((v)/8.0)+15.0*exp((v)/50.0));
        tauqa = 1.0/(alphaqa + betaqa);
        alphaqi = 0.15*1.0/(3100.0*exp((v+10.0)/13.0)+700.3*exp((v+10.0)/70.0));
        betaqi =  0.15*1.0/(95.7*exp(-(v+10.0)/10.0) + 50.0*exp(-(v+10.0)/700.0)) + 0.000229/(1+exp(-(v+10.0)/5.0));
        qi = alphaqi/(alphaqi + betaqi);
        tauqi = 1.0/(alphaqi + betaqi);

        double ist = gst*dst*fst*(v - eist);

    /* Ib ************************************************************************/

        double ibca,ibna,ibk;
        ibna = gbna*(v - ena);
        ibca = gbca*(v - eca);
        ibk  = gbk*(v - ek);

        double ib = (ibna + ibca + ibk);

    /*IK1**********************************************************************/

        double xk1inf;
        xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));

        double ik1 = gk1*xk1inf*(ko/(ko + 0.228880))*(v - ek);

    /**ICaT Cav3.1**************************************************************/

        double dt_inf, tau_dt, ft_inf, tau_ft;
        tau_dt = 1.0/(1.068*exp((v + 26.3)/30.0) + 1.068*exp(-(v + 26.3)/30.0));
        dt_inf = 1.0/(1.0+exp(-(v + 26.0)/6.0));
        tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38));

        ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6));

        double icat = gcat*ft*dt*(v - ecat);

    ///*Ikr********************************************************************/
    ////////////////Changed by wei 23092020/////////
        double ikr_act_inf_f, tau_ikr_act_f, ikr_inact_inf, tau_ikr_inact, ikr_act_inf_s, tau_ikr_act_s, FFF;
    //////////////////////////////////////////////TAU_f//////////////////////////////////////
        ikr_act_inf_f = 1.0/(1.0 + exp(-(v+21.173694)/9.757086));
    //	  p0= [ 1.36363765e+01  7.36947104e-03 -1.48435756e+01  1.29312340e-01
    //	   -9.49390754e+00]
        tau_ikr_act_f = 1.36363765e+01/( 7.36947104e-03*exp((v)/-1.48435756e+01) +  1.29312340e-01*exp(-(v)/-9.49390754e+00));
    //////////////////////////////////////////////TAU_S//////////////////////////////////////
        ikr_act_inf_s = 1.0/(1.0 + exp(-(v+21.173694)/9.757086));
        tau_ikr_act_s = 3.96851039e+02/(5.11087460e+00*exp((v)/1.13591973e+01) + 8.24877522e-02*exp(-(v)/2.37122583e+01));
    //////////////////////////////////////////////Inactivate//////////////////////////////////////
        ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0)));
        tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));
    ////p0= [  6.24211258   4.8009683  -50.31969009  57.62015753 -14.47035665]
        FFF =   6.24211258 /(4.8009683 *exp((v)/(-50.31969009)) +57.62015753 *exp(-(v)/(-14.47035665)));
        ikr_act = (ikr_act_f*(1/(FFF+1)))+(ikr_act_s*(FFF/(FFF+1)));

        double ikr = (1* gkr)*ikr_act*ikr_inact*(v - ek);

    /**IKs********************************************************************/

        double iks_act_inf,tau_iks_act;
        iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
        tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539));

        double iks = gks*iks_act*iks_act*(v - eks);

    /*ICaL*******************************************************************/

        double alpha_dl, beta_dl, tau_dl, tau_fl12, tau_fl13;
        double dl12_inf, fl12_inf;
        double dl13_inf, fl13_inf;
        double fca_inf, taufca;
        if(fabs(v)<=0.001)
        {
            alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)+408.173;
        }
        else if(fabs(v+35.0)<=0.001)
        {
            alpha_dl  = 70.975-84.9*v/(exp(-0.208*v)-1.0);
        }
        else if(fabs(v)>0.001&&fabs(v+35.0)>0.001)
        {
            alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)-84.9*v/(exp(-0.208*v)-1.0);
        }

        if(fabs(v-5.0)<=0.001)
            beta_dl   = 28.575;
        else if(fabs(v-5.0)>0.001)
            beta_dl   = 11.43*(v-5.0)/(exp(0.4*(v-5.0))-1.0);

        tau_dl  = 2000.0/(alpha_dl +beta_dl);
        dl13_inf = 1.0/(1+exp(-(v+13.5)/6.0));
        fl13_inf = 1.0/(1+exp((v+35.0)/7.3));
        tau_fl13 = (7.4 + 45.77*exp(-0.5*(v+28.7)*(v+28.7)/(11*11)));
        dl12_inf = 1.0/(1+exp(-(v+3.0)/5.0));
        fl12_inf = 1.0/(1+exp((v+36.0)/4.6));
        tau_fl12 = (7.4 + 45.77*exp(-0.5*(v+24.7)*(v+24.7)/(11*11)));
        fca_inf = kmfca/(kmfca+casub);
        taufca = fca_inf/alpha_fca;


        double ical12 = SAN_config.ICaL_Block*gcal12*fl12*dl12*fca*(v-ecal);
        double ical13 = SAN_config.ICaL_Block*gcal13*fl13*dl13*fca*(v-ecal);
    //    ica_new = ical12 + ical13;

    /**INa**********************************************************************/

        double m3_inf_ttxr, h_inf_ttxr;
        double m3_inf_ttxs, h_inf_ttxs;
        double m_inf_ttxr ,m_inf_ttxs;
        double tau_mr, tau_hr, tau_jr;
        double fna, hs, hsr;
        fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2;

        m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));

        h_inf_ttxr = 1.0/(1.0 + exp((v+62.578120 )/6.084036));

    //    ///////////////////h_inf_ttxr by wei///////////////////////
    //
    //    h_inf_ttxr = 1.0/(1.0 + exp((v+62.578120+5.0 )/6.084036));

        m3_inf_ttxs = 1.0/(1.0 + exp(-(v+31.097331)/5.0));
        h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));
        m_inf_ttxr = pow(m3_inf_ttxr,0.333);
        m_inf_ttxs = pow(m3_inf_ttxs,0.333);
        tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
        tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
        tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
        hs = (1.0-fna)*h_ttxs+fna*j_ttxs;
        hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;

        double ina_ttxs;
        double ina_ttxr;
        if(fabs(v)>0.005) {
            ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
        } else {
            ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
        }
        if(fabs(v)>0.005) {
            ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
        } else {
            ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));
        }
        ina_ttxr = SAN_config.INa_Block*ina_ttxr; //Apply channel blocks
        ina_ttxs = SAN_config.INa_Block*ina_ttxs;
        double ina = ina_ttxs + ina_ttxr;

        //current = m_ttxr*m_ttxr*m_ttxr;

    /**If**************************************************************************/

        double ihk, ihna, ihk_hcn1, ihk_hcn2, ihk_hcn4, ihna_hcn1, ihna_hcn2, ihna_hcn4;
        double y_inf_1, y_inf_2, y_inf_4, tau_y_1, tau_y_2, tau_y_4;

        y_inf_1 = 1.0/(1.0 + exp((v+75.2)/6.57));
        y_inf_2 = 1.0/(1.0 + exp((v+92.0)/9.5));   //y_inf_2 = 1.0/(1.0 + exp((v+102.0)/5.49));
        y_inf_4 = 1.0/(1.0 + exp((v+91.2)/10.5));  //y_inf_4 = 1.0/(1.0 + exp((v+101.0)/9.71));
        tau_y_1 = 0.0035017700/(exp(-(v+590.192)*0.0208819)+ exp((v-85.4612)/13.26));
        tau_y_2 = 437.116/(0.000225*exp(-v/13.02)+105.395*exp(v/13.02));
        tau_y_4 = 3121.55/(2.758e-05*exp(-v/9.92)+766.3*exp(v/9.92));

        ihk_hcn1 = ghk1*y_1*(v-ek);
        ihk_hcn2 = ghk2*y_2*(v-ek);
        ihk_hcn4 = ghk4*y_4*(v-ek);
        ihna_hcn1 = ghna1*y_1*(v-ena);
        ihna_hcn2 = ghna2*y_2*(v-ena);
        ihna_hcn4 = ghna4*y_4*(v-ena);

        ihk = ihk_hcn1+ihk_hcn2+ihk_hcn4;
        ihna = ihna_hcn1+ihna_hcn2+ihna_hcn4;

    ////////////////Original ih//////////
        double ih = (ihk + ihna);

        //current=ihk_hcn1+ihna_hcn1;
        //=ihk_hcn2+ihna_hcn2;
        //=ihk_hcn4+ihna_hcn4;
    /*Ito*************************************************************************/

        double r_inf, tau_r, q_inf,tau_q;
        q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
        tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67;
        r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
        tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;

        double ito = gto*q*r*(v-ek);

    /*Isus***********************************************************************/

        double isus = gsus*r*(v-ek);

    /*Inak***********************************************************************/

        double inak = inakmax*(pow(ko,1.2)/(pow(kmkp,1.2)+pow(ko,1.2)))*(pow(nai,1.3)/(pow(kmnap,1.3)+pow(nai,1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));

    /****iNaCa*******************************************************************/

        double di,doo,k43,k12,k14,k41,k34,k21,k23,k32,x1,x2,x3,x4;
        di=1+(casub/Kci)*(1+exp(-Qci*v*F/(R*T))+nai/Kcni)+(nai/K1ni)*(1+(nai/K2ni)*(1+nai/K3ni));
        doo=1+(cao/Kco)*(1+exp(Qco*v*F/(R*T)))+(nao/K1no)*(1+(nao/K2no)*(1+nao/K3no));
        k43=nai/(K3ni+nai);
        k12=(casub/Kci)*exp(-Qci*v*F/(R*T))/di;
        k14=(nai/K1ni)*(nai/K2ni)*(1+nai/K3ni)*exp(Qn*v*F/(R*T*2.0))/di;
        k41=exp(-Qn*v*F/(R*T*2.0));
        k34=nao/(K3no+nao);
        k21=(cao/Kco)*exp(Qco*v*F/(R*T))/doo;
        k23=(nao/K1no)*(nao/K2no)*(1+nao/K3no)*exp(-Qn*v*F/(R*T*2.0))/doo;
        k32=exp(Qn*v*F/(R*T*2));
        x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
        x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
        x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
        x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);

        double inaca = kNaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);

    /****Ca handling in the SR***************************************************/

        double Jrel, Jup, Jtr, kcasr, kosrca, kisrca, drdt, dodt, didt, dridt;
        Jup = Pup*(pow(cai/pumpkmf,pumphill) - pow(caup/pumpkmr,pumphill))/(1.0 + pow(cai/pumpkmf,pumphill) + pow(caup/pumpkmr,pumphill));
        Jtr  = (caup - carel)/Ttr;
        Jrel = ks*open*(carel - casub);
        kcasr = maxsr - (maxsr - minsr)/(1.0 + pow(eca50sr/carel,hsrr));
        kosrca = koca/kcasr;
        kisrca = kica*kcasr;

        drdt = kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open;
        dodt = kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated;
        didt = kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated;
        dridt = kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting;

    /****Ca diffusion ***********************************************************/

        double Jcadif;

        Jcadif = (casub - cai)/tdifca;

    /****Ca buffering ***********************************************************/

        double dFtcdt, dFtmcdt, dFtmmdt, dFcmsdt, dFcmidt, dFcqdt;

        dFtcdt = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
        dFtmcdt = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
        dFtmmdt = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
        dFcmsdt = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
        dFcmidt = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
        dFcqdt = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;

    /****Intracellular ionic concentrations *************************************/

        double ca_flux, dcasubdt, dcaidt, dcareldt, dcaupdt, nai_tot, ki_tot;

        ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);
        dcasubdt = (-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcmsdt;
        dcaidt = (Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmidt + ConcTC*dFtcdt + ConcTMC*dFtmcdt);
        dcareldt = Jtr - Jrel - ConcCQ*dFcqdt;
        dcaupdt = Jup-Jtr*vrel/vup;
        nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
        ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;


    /*****Membrane popential ****************************************************/
        //Outward = iks+ikr+ik1+ibK+inak+isus+ibnac+ito

        double Itot = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito; //[nA]
        double dvdt = (-Itot + I_j)/capacitance; // [V/s]

    /*****ALL Inward current ****************************************************/

        double I_inwardcurrent = ih+ina_ttxr+ina_ttxs+ical12+ical13+ist+icat+inaca+ibna+ibca;
    
    //Obtain rates
    states_dot.v = dvdt;
    states_dot.dst = (qa-dst)/tauqa;
    states_dot.fst = (qi-fst)/tauqi;
    states_dot.dt = (dt_inf - dt)/tau_dt;
    states_dot.ft = (ft_inf - ft)/tau_ft;
    //states_dot.ikr_act = (ikr_act_inf-ikr_act)/tau_ikr_act; 
    states_dot.ikr_act_f = (ikr_act_inf_f-ikr_act_f)/tau_ikr_act_f;
    states_dot.ikr_act_s = (ikr_act_inf_s-ikr_act_s)/tau_ikr_act_s;
    states_dot.ikr_inact = (ikr_inact_inf - ikr_inact)/tau_ikr_inact;
    states_dot.iks_act = (iks_act_inf - iks_act)/tau_iks_act;
    states_dot.fl12 = (fl12_inf - fl12)/tau_fl12;
    states_dot.dl12 = (dl12_inf - dl12)/tau_dl;
    states_dot.fl13 = (fl13_inf - fl13)/tau_fl13;
    states_dot.dl13 = (dl13_inf - dl13)/tau_dl;
    states_dot.r = (r_inf-r)/tau_r;
    states_dot.m_ttxr = (m_inf_ttxr - m_ttxr)/tau_mr;
    states_dot.h_ttxr = (h_inf_ttxr - h_ttxr)/tau_hr;
    states_dot.j_ttxr = (h_inf_ttxr - j_ttxr)/tau_jr;
    states_dot.m_ttxs = (m_inf_ttxs - m_ttxs)/tau_mr;
    states_dot.h_ttxs = (h_inf_ttxs - h_ttxs)/tau_hr;
    states_dot.j_ttxs = (h_inf_ttxs - j_ttxs)/tau_jr;
    states_dot.y_1 = (y_inf_1 - y_1)/tau_y_1;
    states_dot.carel = dcareldt;
    states_dot.caup = dcaupdt;
    states_dot.casub = dcasubdt;
    states_dot.Ftc = dFtcdt;
    states_dot.Ftmc = dFtmcdt;
    states_dot.Ftmm = dFtmmdt;
    states_dot.Fcms = dFcmsdt;
    states_dot.Fcmi = dFcmidt;
    states_dot.Fcq = dFcqdt;
    states_dot.cai = dcaidt;
    states_dot.q = (q_inf-q)/tau_q;
    states_dot.fca = (fca_inf - fca)/taufca;
    states_dot.nai = (nai_tot)/(F*vi);
    states_dot.ki = (ki_tot)/(F*vi);
    states_dot.resting = drdt;
    states_dot.open = dodt;
    states_dot.inactivated = didt;
    states_dot.resting_inactivated = dridt;
    states_dot.y_2 = (y_inf_2 - y_2)/tau_y_2;
    states_dot.y_4 = (y_inf_4 - y_4)/tau_y_4;

    //Update outputs
    dependents.ist = ist;
    dependents.ib = ib;
    dependents.ik1 = ik1;
    dependents.icat = icat;
    dependents.ikr = ikr;
    dependents.iks = iks;
    dependents.ical12 = ical12;
    dependents.ical13 = ical13;
    dependents.ina_ttxs = ina_ttxs;
    dependents.ina_ttxr = ina_ttxr;
    dependents.ih = ih;
    dependents.ito = ito;
    dependents.isus = isus;
    dependents.inak = inak;
    dependents.inaca = inaca;
    dependents.ina = ina;
    dependents.I_inwardcurrent = I_inwardcurrent;
    dependents.Itot = Itot;
    dependents.I_j = I_j;
}