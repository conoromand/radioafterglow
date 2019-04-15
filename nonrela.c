/* Updates from ver.1 */
/* The thermal Lorentz factor of the injection electron spectrum is now parameterized by fac_T. */
/* The effect of inverse Compton emission is included. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "dcbh.h"

/* input parameters */

/* ejecta */
const double M_EJ = 1.0e-4*1.99e33; /* total mass [g] */
const double V_EJ_MIN = 0.1*3.0e10; /* minimum initial velcocity [cm/s] */
const double V_EJ_MAX = 0.2*3.0e10; /* maximum initial velocity [cm/s] */
const double POW_EJ_EK = 3.75; /* power low index of the cumulative kinetic energy distribution  */
const int FLAG_EK_PROF = 1; /* frag for the kinetic energy distribution; 0 = monoenergetic ejecta, 1 = power law */

/* ambient medium */
const double RHO_CSM = 10.0*1.67e-24; /* circumstellar mass density [g/cc] */
const double MDOT = 1.0e-5*1.99e33/3.15e7; /* wind mass loss rate of the progenitor [g/s] */
const double V_W = 1.0e8; /* wind velocity [cm/s] */
const double CS_CMS = 3.0e7; /* the sound velocity in the CMS [cm/s] */
const double GAM = 5.0/3.0; /* the adiabatic index of the gas */
const int FLAG_AMB = 1; /* flag for the ambient density strucuture; -1 = wind, 0 = constant, 1 = mix */

/* magnetic field amplification and shock acceleration */
const double EPS_B = 0.33;
const double ZETA_E = 0.4;
const double FRAC_E = 0.1;
const double XI_T = 0.24;
const double POW_ELE = 2.0;

/* time window */
const double t_MIN = 1.0e2; /* in unit of [s] */
const double t_MAX = 1.0*365.0*24.0*60.0*60.0; /* in unit of [s] */
const double TIME_1 = 0.01*24.0*60*60; /* plot time 1 */
const double TIME_2 = 0.1*24.0*60*60; /* plot time 2 */
const double TIME_3 = 1.0*24.0*60*60; /* plot time 3 */
const double TIME_4 = 10.0*24.0*60*60; /* plot time 4 */
const double TIME_5 = 100.0*24.0*60*60; /* plot time 5 */
const int N_t = 1000; /* divide the time window t_MIN < t < t_MAX with N_t */
const double TIME_RESO = 0.0001; /* resolution of time in terms of dynamical time of the forward shock */

/* electron and photon spectrum */
const int N_ne = 256;
const double GAMMA_ELEC_MIN = 1.0;/* minimum Lorentz factor of electron */
const double GAMMA_ELEC_MAX = 2.0e9;/* maximum Lorentz factor of electron */
const int N_nph = 256; /* need to be the same as N_ne */
const double ENE_PH_MIN_DIV_MeC2 = 1.0e-14;/* mimimum energy of photon in unit of electron mass */
const double ENE_PH_MAX_DIV_MeC2 = 1.0;/* maximum energy of photon in unit of electron mass */
const int N_mu = 32;
const double NU_1 = 1.0e8; /* plot frequency 1 */
const double NU_2 = 1.0e9; /* plot frequency 2 */
const double NU_3 = 1.0e10; /* plot frequency 3 */
const double NU_4 = 1.0e11; /* plot frequency 4 */
const double NU_5 = 1.5e18; /* plot frequency 5 */

void shock_dynamics(double r, double *dt, double *dr, double *dr_sh, double *v, double *rho, double *rho_sh, double *Vol_sh);
void nth_e_and_sh_B( double t, double v, double rho,double *B_sh, double *gamma_e_inj, double *gamma_e_max, double *gamma_e_cool, double *gamma_e_th);
void elec_injection(double r, double v, double rho_sh, double B_sh, double gamma_e_inj, double gamma_e_th, double gamma_e_max, double gamma_e[], double dne_dt_inj[]);
double power_ad(double gamma_e, double r, double v);
double power_syn(double gamma_e, double B_sh);
void elec_cooling(double r, double v, double B_sh, double t_ad[], double t_syn[], double P_ad[], double P_syn[], double P_cool[], double gamma_e[]);
void elec_time_evolution(double dt, double gamma_e[], double dene_e[], double ne_old[], double ne_new[], double dne_dt_inj[], double P_cool[]);
double syn_func_fit(double x);
void syn_spec(double B_sh, double dr_sh, double Vol_sh, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double mu[], double del_mu, double gamma_ph[], 
	      double P_nu_syn[], double alpha_nu_syn[], double tau_sa[], double P_nu_syn_esc[]);
void ic_spec(double Vol_sh, double dr_sh, double P_nu_syn_esc[], double ne[], double gamma_e[], double del_ln_gamma_e, double gamma_ph[], double del_ln_gamma_ph, double P_nu_ic[]);
void set_var(double mu[], double *del_mu, double gamma_e[], double dene_e[], double *del_ln_gamma_e, double gamma_ph[], double dene_ph[], double *del_ln_gamma_ph);

int main()
{
  /* setting variables */
  double mu[N_mu],gamma_e[N_ne],dene_e[N_ne],gamma_ph[N_nph],dene_ph[N_nph];
  double del_mu=0.0,del_ln_gamma_e=0.0,del_ln_gamma_ph=0.0;
  set_var(mu,&del_mu,gamma_e,dene_e,&del_ln_gamma_e,gamma_ph,dene_ph,&del_ln_gamma_ph);
  int i,j;

  double t=0.0,dt=0.0,r=0.0,dr=0.0,dr_sh=0.0,v=0.0,rho=0.0,rho_sh=0.0,Vol_sh=0.0;
  double gamma_e_inj=0.0,B_sh=0.0,gamma_e_max=0.0,gamma_e_cool=0.0,gamma_e_th=0.0;
  double norm_th = 0.0, norm_nth = 0.0;
  double dne_dt_inj[N_ne],t_syn[N_ne],t_ad[N_ne],P_syn[N_ne],P_ad[N_ne],P_cool[N_ne],ne_old[N_ne],ne_new[N_ne];
  double P_nu_syn[N_nph],alpha_nu_syn[N_nph],tau_sa[N_nph],P_nu_syn_esc[N_nph],P_nu_ic[N_nph];

  double E_sh_tot=0.0,N_sh_tot=0.0,E_e_tot=0.0,N_e_tot=0.0,E_cool_tot=0.0;

  /* time evolution */
  double del_ln_t = (log(t_MAX/t_MIN))/(double)(N_t-1);
   
  for (i=0;i<N_ne;i++) {
    ne_old[i] = 0.0;
    ne_new[i] = 0.0;
  }
    
  for (i=0;i<N_nph;i++) {
    P_nu_syn[i] = 0.0;
    alpha_nu_syn[i] = 0.0;
  }
    
    int time1 = (int)(log(TIME_1/t_MIN)/del_ln_t);
    int time2 = (int)(log(TIME_2/t_MIN)/del_ln_t);
    int time3 = (int)(log(TIME_3/t_MIN)/del_ln_t);
    int time4 = (int)(log(TIME_4/t_MIN)/del_ln_t);
    int time5 = (int)(log(TIME_5/t_MIN)/del_ln_t);

    int timei = 0;
    int plot_num_t = 0;
    int flag_time = 1;
    int loop_cnt = 0;

    int nu1 = (int)(log(H*NU_1/MeC2/ENE_PH_MIN_DIV_MeC2)/del_ln_gamma_ph);
    int nu2 = (int)(log(H*NU_2/MeC2/ENE_PH_MIN_DIV_MeC2)/del_ln_gamma_ph);
    int nu3 = (int)(log(H*NU_3/MeC2/ENE_PH_MIN_DIV_MeC2)/del_ln_gamma_ph);
    int nu4 = (int)(log(H*NU_4/MeC2/ENE_PH_MIN_DIV_MeC2)/del_ln_gamma_ph);
    int nu5 = (int)(log(H*NU_5/MeC2/ENE_PH_MIN_DIV_MeC2)/del_ln_gamma_ph);
    
    FILE *op1,*op2,*op3;
    char file_name1[256] = {"\0"};
    char file_name2[256] = {"\0"};
    char file_name3[256] = {"\0"};
    char dat[256] = ".dat";
    
    sprintf(file_name1,"dynamics_cygnus%s",dat);
    op1 = fopen(file_name1,"w+");
    fprintf(op1,"# t [s], r [cm], v_sh[cm/s], rho_sh [g/cm], Vol_sh [cm^3], dr_sh [cm], B_sh [G], gamma_e_inj \n");
    
    sprintf(file_name2,"lightcurve_cygnus%s",dat);
    op2 = fopen(file_name2,"w+");
    fprintf(op2,"# t [s], L_nu @ %12.3e [Hz], %12.3e [Hz], %12.3e [Hz], %12.3e [Hz], %12.3e [Hz] \n",
            gamma_ph[nu1]*MeC2/H,gamma_ph[nu2]*MeC2/H,gamma_ph[nu3]*MeC2/H,gamma_ph[nu4]*MeC2/H,gamma_ph[nu5]*MeC2/H);
    
    t = t_MIN;
    if (FLAG_EK_PROF == 0)
        r = V_EJ_MIN*t;
    else
        r = V_EJ_MAX*t;
    
    while (t < t_MAX) {
        timei = (int)(log(t/t_MIN)/del_ln_t);
    
        shock_dynamics(r,&dt,&dr,&dr_sh,&v,&rho,&rho_sh,&Vol_sh);
        nth_e_and_sh_B(t,v,rho,&B_sh,&gamma_e_inj,&gamma_e_max,&gamma_e_cool,&gamma_e_th);
        elec_injection(r,v,rho,B_sh,gamma_e_inj,gamma_e_th,gamma_e_max,gamma_e,dne_dt_inj);
        elec_cooling(r,v,B_sh,t_ad,t_syn,P_ad,P_syn,P_cool,gamma_e);
        elec_time_evolution(dt,gamma_e,dene_e,ne_old,ne_new,dne_dt_inj,P_cool);
        syn_spec(B_sh,dr_sh,Vol_sh,ne_new,gamma_e,dene_e,del_ln_gamma_e,mu,del_mu,gamma_ph,P_nu_syn,alpha_nu_syn,tau_sa,P_nu_syn_esc);

        for (i=0;i<N_ne;i++){
            ne_old[i] = ne_new[i];
        }
    
        if (loop_cnt % (int)(1.0/(TIME_RESO*10.0)) == 0){
            ic_spec(Vol_sh,dr_sh,P_nu_syn_esc,ne_new,gamma_e,del_ln_gamma_e,gamma_ph,del_ln_gamma_ph,P_nu_ic);
            fprintf(op1,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",t,r,v,rho,Vol_sh,dr_sh,B_sh,gamma_e_inj);
            fprintf(op2,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
                    t,
                    (P_nu_syn_esc[nu1]+P_nu_ic[nu1])*Vol_sh/4.0/M_PI,
                    (P_nu_syn_esc[nu2]+P_nu_ic[nu2])*Vol_sh/4.0/M_PI,
                    (P_nu_syn_esc[nu3]+P_nu_ic[nu3])*Vol_sh/4.0/M_PI,
                    (P_nu_syn_esc[nu4]+P_nu_ic[nu4])*Vol_sh/4.0/M_PI,
                    (P_nu_syn_esc[nu5]+P_nu_ic[nu5])*Vol_sh/4.0/M_PI);
        }
        
        if ((timei == time1 && flag_time == 1) || (timei == time2 && flag_time == 2) || (timei == time3 && flag_time == 3)
            || (timei == time4 && flag_time == 4) || (timei == time5 && flag_time == 5)){
            plot_num_t += 1;
            ic_spec(Vol_sh,dr_sh,P_nu_syn_esc,ne_new,gamma_e,del_ln_gamma_e,gamma_ph,del_ln_gamma_ph,P_nu_ic);
            sprintf(file_name3,"Lnu_cygnus%d%s",plot_num_t,dat);
            op3 = fopen(file_name3,"w+");
            fprintf(op3,"# t = %12.3e [day] \n",t/(24.0*60.0*60.0));
            fprintf(op3,"# nu [Hz], P_nu [erg/s/cm^3/Hz], tau_nu, P_nu_esc [erg/s/cm^3/Hz], L_nu_syn [erg/s/Hz], L_nu_ic [erg/s/Hz], L_nu_tot [erg/s/Hz], \n");
            for (j=0;j<N_nph;j++) {
                fprintf(op3,"%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e \n",
                        gamma_ph[j]*MeC2/H,P_nu_syn[j],tau_sa[j],P_nu_syn_esc[j],P_nu_syn_esc[j]*Vol_sh/4.0/M_PI,P_nu_ic[j]*Vol_sh/4.0/M_PI,(P_nu_syn_esc[j]+P_nu_ic[j])*Vol_sh/4.0/M_PI);
            }
            fclose(op3);
            
            sprintf(file_name3,"ne_cygnus%d%s",plot_num_t,dat);
            op3 = fopen(file_name3,"w+");
            fprintf(op3,"# t = %12.3e [day] \n",t/(24.0*60.0*60.0));
            fprintf(op3,"# Gamma_e, n_e [1/cm^3] \n");
            for (j=0;j<N_nph;j++) {
                fprintf(op3,"%12.3e %12.3e %12.3e \n",gamma_e[j],gamma_e[j]*ne_new[j]*dene_e[j],ne_new[j]/Vol_sh*dene_e[j]/MeC2);
            }
            fclose(op3);
            flag_time += 1;
        }

        /* checking the energy and number conservation */
        N_sh_tot += rho*(4.0*M_PI*r*r*dr)/M_PRO;
        E_sh_tot += 0.5*gamma_e[0]*MeC2*dne_dt_inj[0]*dt*gamma_e[0]*del_ln_gamma_e;
        N_e_tot = 0.5*ne_new[0]*gamma_e[0]*del_ln_gamma_e;
        E_e_tot = 0.5*gamma_e[0]*MeC2*ne_new[0]*gamma_e[0]*del_ln_gamma_e;
        E_cool_tot += 0.5*P_cool[0]*dt*ne_new[0]*gamma_e[0]*del_ln_gamma_e;
        for (j=1;j<N_ne-1;j++){
            E_sh_tot += gamma_e[j]*MeC2*dne_dt_inj[j]*dt*gamma_e[j]*del_ln_gamma_e;
            N_e_tot += ne_new[j]*gamma_e[j]*del_ln_gamma_e;
            E_e_tot += gamma_e[j]*MeC2*ne_new[j]*gamma_e[j]*del_ln_gamma_e;
            E_cool_tot += P_cool[j]*dt*ne_new[j]*gamma_e[j]*del_ln_gamma_e;
        }
        E_sh_tot += 0.5*gamma_e[N_ne-1]*MeC2*dne_dt_inj[N_ne-1]*dt*gamma_e[N_ne-1]*del_ln_gamma_e;
        N_e_tot += 0.5*ne_new[N_ne-1]*gamma_e[N_ne-1]*del_ln_gamma_e;
        E_e_tot += 0.5*gamma_e[N_ne-1]*MeC2*ne_new[N_ne-1]*gamma_e[N_ne-1]*del_ln_gamma_e;
        E_cool_tot += 0.5*P_cool[N_ne-1]*dt*ne_new[N_ne-1]*gamma_e[N_ne-1]*del_ln_gamma_e;

        printf("t = %12.3e [s] \n",t);
        printf("Number conservation: %12.3e \n",N_e_tot/N_sh_tot);
        printf("Energy conservation: %12.3e \n",(E_e_tot+E_cool_tot)/E_sh_tot);
        printf("\n");
        
        t += dt;
        r += dr;
        loop_cnt += 1;

    }
    fclose(op1);
    fclose(op2);

    return 0;
}

void shock_dynamics(double r, double *dt, double *dr, double *dr_sh, double *v, double *rho, double *rho_sh, double *Vol_sh)
{
  double dt_tmp=0.0,dr_tmp=0.0,dr_sh_tmp=0.0,rho_tmp=0.0,rho_sh_tmp=0.0,v_tmp=0.0,Vol_sh_tmp=0.0;
  double M_sh_tot;
  double v_refreshed_sh = 0.0;

    /* assuming a compression ratio of 4 -> strong shock condition */
    if (FLAG_AMB == -1){
        rho_sh_tmp = MDOT/M_PI/V_W/r/r;
        M_sh_tot = MDOT/V_W*r;
    } else if (FLAG_AMB == 0){
        rho_sh_tmp = 4.0*RHO_CSM;
        M_sh_tot = RHO_CSM*4.0*M_PI/3.0*r*r*r;
    } else if (FLAG_AMB == 1){
        double r_w = sqrt(GAM*MDOT*V_W/8.0/M_PI/RHO_CSM/CS_CMS/CS_CMS);
        if (r < r_w){
            rho_sh_tmp = MDOT/M_PI/V_W/r/r;
            M_sh_tot = MDOT/V_W*r;
        } else {
            rho_sh_tmp = 4.0*RHO_CSM;
            M_sh_tot = MDOT/V_W*r_w+RHO_CSM*4.0*M_PI/3.0*(r*r*r-r_w*r_w*r_w);
        }
    }
    rho_tmp = rho_sh_tmp/4.0; 

    if (FLAG_EK_PROF == 0)
      v_tmp = sqrt(M_EJ/(M_EJ+M_sh_tot))*V_EJ_MIN;
    else {
      v_refreshed_sh = pow(M_EJ/M_sh_tot,1.0/(POW_EJ_EK+2.0))*V_EJ_MIN;
      if (v_refreshed_sh > V_EJ_MAX)
	v_tmp = V_EJ_MAX;
      else { 
	if (M_EJ > M_sh_tot)
            v_tmp = pow(M_EJ/M_sh_tot,1.0/(POW_EJ_EK+2.0))*V_EJ_MIN;
        else
            v_tmp = sqrt(M_EJ/M_sh_tot)*V_EJ_MIN;
      }
    }
    
    Vol_sh_tmp = M_sh_tot/rho_sh_tmp;
    dr_sh_tmp = Vol_sh_tmp/(4.0*M_PI*r*r);
    dt_tmp = TIME_RESO*r/v_tmp;
    dr_tmp = v_tmp*dt_tmp;
    
    *dt = dt_tmp;
    *dr = dr_tmp;
    *dr_sh = dr_sh_tmp;
    *v = v_tmp;
    *rho = rho_tmp;
    *rho_sh = rho_sh_tmp;
    *Vol_sh = Vol_sh_tmp;
}

void nth_e_and_sh_B(double t, double v, double rho, double *B_sh, double *gamma_e_inj, double *gamma_e_max, double *gamma_e_cool, double *gamma_e_th)
{
    double B_sh_tmp = sqrt(9.0*M_PI*EPS_B*rho*v*v);

    *B_sh = B_sh_tmp;
    *gamma_e_inj = 0.5*ZETA_E*M_PRO/M_ELE*pow(v/C,2.0);
    *gamma_e_max = sqrt(9.0*M_PI*ELEC*v*v/10.0/SIGMA_T/B_sh_tmp/C/C);
    *gamma_e_cool = 6.0*M_PI*M_ELE*C/SIGMA_T/B_sh_tmp/B_sh_tmp/t;
    *gamma_e_th = 0.5*XI_T*M_PRO/M_ELE*pow(v/C,2.0);
}


void elec_injection(double r, double v, double rho, double B_sh, double gamma_e_inj, double gamma_e_th, double gamma_e_max, double gamma_e[], double dne_dt_inj[])
{
  double integ_th=0.0,integ_nth=0.0,gamma_e_th_tmp=1.0;
  double del_ln_gamma_th = log(gamma_e_max)/(double)(N_ne-1);
  int i;

  for (i=1;i<N_ne-1;i++)
    {
      gamma_e_th_tmp = exp(del_ln_gamma_th*(double)i);
      integ_th += (gamma_e_th_tmp*sqrt(gamma_e_th_tmp*gamma_e_th_tmp-1.0)*exp(-gamma_e_th_tmp/gamma_e_th)*gamma_e_th_tmp)*del_ln_gamma_th;
    }
  integ_th += 0.5*(gamma_e_max*sqrt(gamma_e_max*gamma_e_max-1.0)*exp(-gamma_e_max/gamma_e_th)*gamma_e_max)*del_ln_gamma_th;

  if (POW_ELE != 1.0){
    integ_nth = (pow(gamma_e_inj,1.0-POW_ELE)-pow(gamma_e_max,1.0-POW_ELE))/(POW_ELE-1.0);
  } else {
    integ_nth = log(gamma_e_max/gamma_e_inj);
  }

  double dN_ele_tot = 4.0*M_PI*r*r*v*rho/M_PRO;
  double norm_th = (1.0-FRAC_E)*dN_ele_tot/integ_th;
  double norm_nth = FRAC_E*dN_ele_tot/integ_nth;

  for (i=0;i<N_ne;i++){
    if (gamma_e[i] > gamma_e_inj && gamma_e[i] < gamma_e_max)
      dne_dt_inj[i] = norm_th*gamma_e[i]*sqrt(gamma_e[i]*gamma_e[i]-1.0)*exp(-gamma_e[i]/gamma_e_th)+norm_nth*pow(gamma_e[i],-POW_ELE);
    else
      dne_dt_inj[i] = norm_th*gamma_e[i]*sqrt(gamma_e[i]*gamma_e[i]-1.0)*exp(-gamma_e[i]/gamma_e_th);
  }
}


double power_ad(double gamma_e, double r, double v)
{
  return gamma_e*MeC2*v/r;
}

double power_syn(double gamma_e, double B_sh)
{
  // electron synchrotron energy loss rate (see e.g., Eq. 7.13 of Dermer & Menon)
  // double sin2phi = 2.0/3.0; /* averaging pitch angle */
  // double beta_par = 1.0; /* assuming that particles are relativistic */
  
  return 4.0/3.0*C*SIGMA_T*(B_sh*B_sh/8.0/M_PI)*gamma_e*gamma_e;
}


void elec_cooling(double r, double v, double B_sh,double t_ad[], double t_syn[], double P_ad[], double P_syn[], double P_cool[], double gamma_e[])
{
  int i;
  for (i=0;i<N_ne;i++) {
    P_ad[i] = power_ad(gamma_e[i],r,v);
    P_syn[i] = power_syn(gamma_e[i],B_sh);
    P_cool[i] = P_ad[i]+P_syn[i];
    t_ad[i] = gamma_e[i]*MeC2/P_ad[i];
    t_syn[i] = gamma_e[i]*MeC2/P_syn[i];
  }
}

void elec_time_evolution(double dt, double gamma_e[], double dene_e[], double ne_old[], double ne_new[], double dne_dt_inj[], double P_cool[])
{
  int i;

  ne_new[N_ne-1] = (ne_old[N_ne-1]+dne_dt_inj[N_ne-1]*dt)/(1.0+dt/dene_e[N_ne-1]*P_cool[N_ne-1]);
  for(i=N_ne-2;i>0;i--){
    ne_new[i] = (ne_old[i]+dne_dt_inj[i]*dt+ne_old[i+1]*dt/dene_e[i]*P_cool[i+1])/(1.0+dt/dene_e[i]*P_cool[i]);
  }
  ne_new[0] = ne_old[0]+dne_dt_inj[0]+ne_old[1]*dt/dene_e[1]*P_cool[1]/(1.0+dt/dene_e[0]*P_cool[0]);

}


/* Synchrotron emission */
double syn_func_fit(double x)
{
    /* analytical fitting of synchrotron function F(x) */
    /* see http://arxiv.org/pdf/1301.6908.pdf */
    
    double F1 = M_PI*pow(2.0,5.0/3.0)/sqrt(3.0)/GAMMA13*pow(x,1.0/3.0);
    double F2 = sqrt(M_PI/2.0)*exp(-x)*pow(x,1.0/2.0);
    
    double a1_1 = -0.97947838884478688;
    double a1_2 = -0.83333239129525072;
    double a1_3 = 0.1554179602681624;
    double H_1 = a1_1*pow(x,1.0)+a1_2*pow(x,1.0/2.0)+a1_3*pow(x,1.0/3.0);
    double delta_1 = exp(H_1);
    
    double a2_1 = -0.0469247165562628882;
    double a2_2 = -0.70055018056462881;
    double a2_3 = 0.0103876297841949544;
    double H_2 = a2_1*pow(x,1.0)+a2_2*pow(x,1.0/2.0)+a2_3*pow(x,1.0/3.0);
    double delta_2 = 1.0-exp(H_2);
    
    return F1*delta_1+F2*delta_2;
}


void syn_spec(double B_sh, double dr_sh, double Vol_sh, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double mu[], double del_mu, double gamma_ph[], 
	      double P_nu_syn[], double alpha_nu_syn[], double tau_sa[], double P_nu_syn_esc[])
{
    int i,j,k;
    double nu,x,sin_alpha;
    double integ=0.0;
    double integ_alpha=0.0;
    
    for (k=0;k<N_nph;k++) {
        nu = gamma_ph[k]*MeC2/H;
        for (i=0;i<N_ne;i++) {
        //for (j=1;j<N_mu-1;j++) {
	      /* Do not take summation for j = 0 and j = N_mu where x = ∞ and syn_func = 0 */
	      //sin_alpha = sqrt(1.0-mu[j]*mu[j]);
	      sin_alpha = 2.0/3.0;
	      x= (2.0*M_PI*nu)/(3.0*ELEC*gamma_e[i]*gamma_e[i]*B_sh/2.0/M_ELE/C*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */
	      if (i==0 || i==N_ne-1) {
		//integ += 0.5*sin_alpha*ne[i]*syn_func_fit(x)*del_mu/2.0*gamma_e[i]*del_ln_gamma_e;
		//integ_alpha += -0.5*sin_alpha*pow(gamma_e[i],2.0)*(-ne[i]/pow(gamma_e[i],2.0))/dene_e[i]*syn_func_fit(x)*del_mu/2.0*gamma_e[i]*del_ln_gamma_e/MeC2;
		integ += 0.5*sin_alpha*ne[i]*syn_func_fit(x)*gamma_e[i]*del_ln_gamma_e;
		integ_alpha += -0.5*sin_alpha*pow(gamma_e[i],2.0)*(-ne[i]/pow(gamma_e[i],2.0))/dene_e[i]*syn_func_fit(x)*gamma_e[i]*del_ln_gamma_e/MeC2;
	      } else {
		//integ += sin_alpha*ne[i]*syn_func_fit(x)*del_mu/2.0*gamma_e[i]*del_ln_gamma_e;
		//integ_alpha += -sin_alpha*pow(gamma_e[i],2.0)*(ne[i+1]/pow(gamma_e[i+1],2.0)-ne[i]/pow(gamma_e[i],2.0))/dene_e[i]*syn_func_fit(x)*del_mu/2.0*gamma_e[i]*del_ln_gamma_e/MeC2;
		integ += sin_alpha*ne[i]*syn_func_fit(x)*gamma_e[i]*del_ln_gamma_e;
		integ_alpha += -sin_alpha*pow(gamma_e[i],2.0)*(ne[i+1]/pow(gamma_e[i+1],2.0)-ne[i]/pow(gamma_e[i],2.0))/dene_e[i]*syn_func_fit(x)*gamma_e[i]*del_ln_gamma_e/MeC2;
	      }
	      //}
        }
        P_nu_syn[k] = sqrt(3.0)*pow(ELEC,3.0)*B_sh/MeC2*integ/Vol_sh; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
        alpha_nu_syn[k] = C*C/8.0/M_PI/nu/nu*sqrt(3.0)*pow(ELEC,3.0)*B_sh*integ_alpha/Vol_sh; /* Eq. (6.52) of Rybicki & Lightman */
        tau_sa[k] = alpha_nu_syn[k]*dr_sh;
        
        if (tau_sa[k] > 1.0e-6)
            P_nu_syn_esc[k] = (1.0-exp(-tau_sa[k]))*P_nu_syn[k]/tau_sa[k];
        else
            P_nu_syn_esc[k] = P_nu_syn[k];
        
        integ = 0.0;
        integ_alpha = 0.0;
    }
}

void ic_spec(double Vol_sh, double dr_sh, double P_nu_syn_esc[], double ne[], double gamma_e[], double del_ln_gamma_e, double gamma_ph[], double del_ln_gamma_ph, double P_nu_ic[])
{
    /* inverse compton scattering of synchrotron photons by electrons "ne[]" at a energy "gamma_ph1" [1/erg/cm^3/s] */
    /* assuming that both CMB and ne[] are isotropic and using the head-on-collision approximation */
    int i,j,k;
    double E1,Gamma_eps,q,F_c,eph;
    double gamma_e_min,gamma_ph_max,gamma_ph_min;
    double integ=0.0;
    double C3H3 = C*C*C*H*H*H;
    
    for (k=0;k<N_nph;k++) {
        gamma_e_min = gamma_ph[k];
        for (i=0;i<N_ne;i++) {
            gamma_ph_max = gamma_ph[k]/(1.0-gamma_ph[k]/gamma_e[i]);
            gamma_ph_min = gamma_ph[k]/4.0/gamma_e[i]/(gamma_e[i]-gamma_ph[k]);
            E1 = gamma_ph[k]/gamma_e[i];
            for (j=0;j<N_nph;j++) {
                
                Gamma_eps = 4.0*gamma_ph[j]*gamma_e[i];
                q = E1/(1.0-E1)/Gamma_eps;
                F_c = (2.0*q*log(q)+(1.0+2.0*q)*(1.0-q)+0.5*pow(Gamma_eps*q,2.0)/(1.0+Gamma_eps*q)*(1.0-q));
                eph = gamma_ph[j]*MeC2;
                /* see Eq. (2.48) of Blumenthal & Gould Rev. Mod. Phys., 42, 237 */
                
                if (gamma_e[i] <= gamma_e_min || gamma_ph[j] <= gamma_ph_min || gamma_ph[j] >= gamma_ph_max) {
                    integ += 0.0;
                } else if (i==0 || i==N_ne-1 || j==0 || j==N_nph-1) {
                    integ += 0.5*ne[i]/gamma_e[i]*(8.0*M_PI*pow(eph,2.0)/(C3H3)/(exp(eph/K_B/T_CMB)-1.0)+(P_nu_syn_esc[j]*dr_sh/C/H/eph))*F_c*del_ln_gamma_ph*del_ln_gamma_e;
                } else {
                    integ += ne[i]/gamma_e[i]*(8.0*M_PI*pow(eph,2.0)/(C3H3)/(exp(eph/K_B/T_CMB)-1.0)+(P_nu_syn_esc[j]*dr_sh/C/H/eph))*F_c*del_ln_gamma_ph*del_ln_gamma_e;
                }
            }
        }
        P_nu_ic[k] = 3.0/4.0*SIGMA_T*C*gamma_ph[k]*MeC2*H*integ/Vol_sh;
        integ = 0.0;
    }
    
}



void set_var(double mu[], double *del_mu, double gamma_e[], double dene_e[], double *del_ln_gamma_e, double gamma_ph[], double dene_ph[], double *del_ln_gamma_ph)
{
    /* setting basic variables as strings */
    
    int i;
    
    /* electron energy in unit of MeC2 */
    double del_ln_gamma_e_tmp = (log(GAMMA_ELEC_MAX)-log(GAMMA_ELEC_MIN))/(double)(N_ne-1);
    for (i=0;i<N_ne;i++){
        gamma_e[i] = GAMMA_ELEC_MIN*exp(del_ln_gamma_e_tmp*(double)i);
        dene_e[i] = gamma_e[i]*MeC2*(exp(del_ln_gamma_e_tmp)-1.0);
    }
    *del_ln_gamma_e = del_ln_gamma_e_tmp;
    
    /* photon energy in unit of MeC2 */
    double del_ln_gamma_ph_tmp = (log(ENE_PH_MAX_DIV_MeC2)-log(ENE_PH_MIN_DIV_MeC2))/(double)(N_nph-1);
    for (i=0;i<N_nph;i++){
        gamma_ph[i] = ENE_PH_MIN_DIV_MeC2*exp(del_ln_gamma_ph_tmp*(double)i);
        dene_ph[i] = gamma_ph[i]*MeC2*(exp(del_ln_gamma_ph_tmp)-1.0);
    }
    *del_ln_gamma_ph = del_ln_gamma_ph_tmp;
    
    /* scattering angle */
    double del_mu_tmp = 2.0/(double)(N_mu-1);
    for (i=0;i<N_mu;i++){
        mu[i] = -1.0+del_mu_tmp*(double)i;
    }
    *del_mu = del_mu_tmp;
    
}
