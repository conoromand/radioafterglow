#ifndef FUNCLIST_H
#define FUNCLIST_H


void calc_conical_model();
void calc_jet();
void calc_shocked_jet();
void calc_sync_map();
void calc_lightcurve(double nuobs, double thetaobs);

void egg_shape(double nuobs, double tobs, double thetaobs, double mu_integ[], double dmu_integ[], double beam_fac[], double vol_fac[], int *time_index_min, int *time_index_max, int nu_index[]);

void set_integ_base_ne(double ne[], double gam_e[], double dene_e[], double *del_ln_gam_e);
void set_integ_base_pnu(double gam_ph[], double dene_ph[], double *del_ln_gam_ph);

void elec_injection(double t_sh, double n_sh, double B_sh, double gam_e_inj, double gam_e_th, double gam_e_max, double gam_e[], double dne_dt_inj[]);
double power_ad(double gam_e, double t);
double power_syn(double gam_e, double B_sh);
void elec_cooling(double t, double B_sh, double P_cool[], double gam_e[]);
void elec_time_evolution(double dt, double gam_e[], double dene_e[], double ne_old[], double ne_new[], double dne_dt_inj[], double P_cool[]);
double syn_func_fit(double x);
void syn_spec(double B, double dr, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double gamma_ph[], double P_nu_syn[], double alpha_nu_syn[]);

void distances(double z, double *d_h, double *d_c, double *d_a, double *d_l);


#endif
