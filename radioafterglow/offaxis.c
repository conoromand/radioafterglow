#include <stdio.h>
#include <math.h>

const double E_c = 4.3e54; /* [erg] */
const double M_c = 22.0*2e33; /* [g] */
const double R_c = 2.0e13; /* [cm] */
const double kappa = 0.35; /* [cm^2 g^{-1}] */
const double T_ion = 6000; /* [K] */

const double C = 2.9979e10; /* [cm] */
const double SIGMA_T = 6.6524e-25; /* [cm^2] */
const double H = 6.6261e-27; /* [erg s] */
const double K_B = 1.3806e-16; /* [erg/K] */
const double ELEC = 4.80320427e-10; /* [cgs] */
const double M_PRO = 1.6726219e-24; /* [g] */
const double M_ELE = 9.10938356e-28; /* [g] */
const double MeC2 = 8.18709e-7; /* electron mass [erg] */

const double H_0 = 67.74*1.0e5/1.0e6; /* [cm/s/pc] */ 
const double OMEGA_M = 0.3089;
const double OMEGA_LAMBDA = 0.6911;
const double OMEGA_IGM = 0.04374;

const double M_SUN = 1.989e33; /* [g]  */
const double PC = 3.086e18; /* [cm] */
const double YR = 365.0*24.0*60.0*60.0; /* [s]  */

const double GAMMA13 = 2.67893; /* Gamma(1/3) */


/////////////////////////////////////////
/* input parameters */
/////////////////////////////////////////
/* Ambient density profile */
const double rfs_min = 1.0e12; /* minimum radius of the free expansion phase in unit of [cm] */
const double v_w_str=2.0e8;
const double Mdot_str=1.0e-4*M_SUN/YR;
const double k=2.0;
const double A=Mdot_str/4.0/M_PI/v_w_str;
//const double k=0.0;
//const double A=0.1*M_PRO;
const double a=1.0; /* see Granot & Piran 2011 */
const double b=0.45; /* see Granot & Piran 2011 */
/* redshift */
const double z = 0.1;
/* jet parameters  */
const double theta_j_0=2.0*M_PI*5.0/360.0;
const double gam_j_0=200.0;
/* shock parameters */
const double eps_B = 0.001;
const double eps_e = 0.01;
const double p = 2.2;
const double ZETA_E = 0.1;//0.4;
const double FRAC_E = 1.0;//0.01;
const double XI_T = 0.24;
const double POW_ELE = 2.0;

const int N_tbin = 500;
const int N_ne = 512;   /* need to be the same as N_nph */
const int N_nph = 512; /* need to be the same as N_ne */
const int N_mu = 64;

const double GAMMA_ELEC_MIN = 1.0;/* minimum Lorentz factor of electron */
const double GAMMA_ELEC_MAX = 2.0e9;/* maximum Lorentz factor of electron */
const double ENE_PH_MIN_DIV_MeC2 = 1.0e-14;/* mimimum energy of photon in unit of electron mass */
const double ENE_PH_MAX_DIV_MeC2 = 1.0e3;/* maximum energy of photon in unit of electron mass */

/* path */
const char path[256]="/Users/kakashi/offaxis_data/";
/////////////////////////////////////////


void calc_conical_model();
void calc_jet();
void calc_shocked_jet();
void calc_sync_map();

void egg_shape(double nuobs, double tobs, double z);

void set_integ_base_ne(double ne[], double gam_e[], double dene_e[], double *del_ln_gam_e);
void set_integ_base_pnu(double gam_ph[], double dene_ph[], double *del_ln_gam_ph);
void set_integ_base_mu(double mu[], double *del_mu);

void elec_injection(double t_sh, double n_sh, double B_sh, double gam_e_inj, double gam_e_th, double gam_e_max, double gam_e[], double dne_dt_inj[]);
double power_ad(double gam_e, double t);
double power_syn(double gam_e, double B_sh);
void elec_cooling(double t, double B_sh, double P_cool[], double gam_e[]);
void elec_time_evolution(double dt, double gam_e[], double dene_e[], double ne_old[], double ne_new[], double dne_dt_inj[], double P_cool[]);
double syn_func_fit(double x);
void syn_spec(double B, double dr, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double gamma_ph[], double P_nu_syn[], double alpha_nu_syn[]);

void distances(double z, double *d_h, double *d_c, double *d_a, double *d_l);


int main()
{
    /*
    calc_conical_model();
    calc_jet();
    calc_shocked_jet();
    calc_sync_map();
    */
     
    /*
    int i,j;
    double ne_list[N_tbin][N_ne],Pnu_list[N_tbin][N_nph];
    
    FILE *ip1,*ip2;
    char head1[256]="map_ne",head2[256]="map_pnu",dat[256]=".dat",input_file_name1[256]={"\0"},input_file_name2[256]={"\0"};
    sprintf(input_file_name1,"%s%s%s",path,head1,dat);
    sprintf(input_file_name2,"%s%s%s",path,head2,dat);
    ip1 = fopen(input_file_name1,"r");
    ip2 = fopen(input_file_name2,"r");
    
    for (i=0; i<N_tbin-1; i++) {
        for (j=0; j<N_ne; j++) {
            fscanf(ip1,"%le ",&ne_list[i][j]);
        }
        fscanf(ip1,"\n");
    }
    fclose(ip1);
    for (i=0; i<N_tbin-1; i++) {
        for (j=0; j<N_nph; j++) {
            fscanf(ip2,"%le ",&Pnu_list[i][j]);
        }
        fscanf(ip2,"\n");
    }
    fclose(ip2);
    */
    
    double nuobs=1.0e9,tobs = 2.0e7,z=.1;
    egg_shape(nuobs,tobs,z);
    
    return 0;
}


//For a given (Tobs,z) and jet dynamics, calculate an "egg shape" retion that contributes to the observed flux
void egg_shape(double nuobs, double tobs, double z)
{
    /* reading the input file of the jet dynamics */
    double t[2*N_tbin],t_s[2*N_tbin],gam_j[2*N_tbin],R[2*N_tbin],R_para[2*N_tbin];
    FILE *ip;
    char head_ip[256]="conical_jet",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    int i=0;
    while (fscanf(ip,"%le %le %le %le %*le %le %*le %*le \n",&t[i],&t_s[i],&gam_j[i],&R[i],&R_para[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    
    int j;
    double mu[N_mu],Regg[N_mu],beam_fac[N_mu],nu_debeam[N_mu];
    int time_index[N_mu],nu_index[N_mu];
    double del_mu;
    set_integ_base_mu(mu,&del_mu);
    for (j=0; j<N_mu; j++) {
        i=0;
        while ((t[i]-mu[j]*R[i]/C) < tobs/(1.+z) && i<2*N_tbin){
            i++;
        }
        //printf("%12.3e \n",(t[i]-mu[j]*R[i]/C)/(tobs/(1.+z)));
        if (i == 0 || i == 2*N_tbin){
            Regg[j] = 0.0;
        } else {
            Regg[j] = R[i];
        }
        time_index[j] = i;
        beam_fac[j] = gam_j[i]*gam_j[i]*pow(1.-sqrt(1.-1./gam_j[i]/gam_j[i])*mu[j],2.0);
        nu_debeam[j] = nuobs*gam_j[i]*(1.-sqrt(1.-1./gam_j[i]/gam_j[i])*mu[j]);
    }
    
    
    FILE *op;
    char head_op[256]="egg",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head_op,dat);
    op = fopen(output_file_name,"w");
    for (i=0; i<N_mu; i++) {
        fprintf(op,"%12.3e %12.3e %d %12.3e %12.3e \n",mu[i],Regg[i],time_index[i],beam_fac[i],nu_debeam[i]);
    }
    fclose(op);
    
}


void calc_conical_model()
{
    FILE *op;
    char head[256]="conical_model",dat[256]=".dat",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head,dat);
    
    int i;
    double r_s_min=pow((1.0-cos(theta_j_0))*(gam_j_0*gam_j_0-1.0),-1.0/(3.0-k)),r_s_max=1.0e1,del_ln_r_s=(log(r_s_max)-log(r_s_min))/(double)(N_tbin-1),r_s=r_s_min,dr_s=0.0;
    double theta_j_tmp=theta_j_0,u_tmp=sqrt(gam_j_0*gam_j_0-1.0);
    double t=r_s_min,dt=0.,t_s=t/(sqrt(u_tmp*u_tmp+1.));
    
    op = fopen(output_file_name,"w+");
    fprintf(op,"#k=%le A=%le a=%le b=%le\n",k,A,a,b);
    for (i=0;i<N_tbin;i++){
        r_s = r_s_min*exp(del_ln_r_s*(double)i);
        dr_s = r_s*(exp(del_ln_r_s)-1.0);
        dt = dr_s*sqrt(u_tmp*u_tmp+1.)/u_tmp;
        t += dt; /* source rest frame */
        t_s += dt/(sqrt(u_tmp*u_tmp+1.)); /* jet rest frame */
        
        theta_j_tmp += 1.0/r_s/pow(1.0+pow(r_s,k-3.0)/(1.0-cos(theta_j_tmp)),(1.0+a)/2.0)/pow(theta_j_tmp,a)*dr_s;
        if (theta_j_tmp>M_PI/2.0)
            theta_j_tmp = M_PI/2.0;
        u_tmp = pow(r_s,-(3.0-k)/2.0)/sqrt(1.0-cos(theta_j_tmp));
        
        fprintf(op,"%le %le %le %le %le %le \n",
                r_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),t*pow(4.0/(3.0-k),-1.0/(3.0-k)),t_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),theta_j_tmp,u_tmp,sqrt(u_tmp*u_tmp+1.0));
    }
    fclose(op);
}


void calc_jet()
{
    /* reading the input file of the jet dynamics */
    double r_nd[N_tbin],t_nd[N_tbin],t_s_nd[N_tbin],theta_j_dc[N_tbin],gam_j_dc[N_tbin];
    FILE *ip;
    char head_ip[256]="conical_model",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    int i=0;
    while (fscanf(ip,"%le %le %le %le %*le %le \n",&r_nd[i],&t_nd[i],&t_s_nd[i],&theta_j_dc[i],&gam_j_dc[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    
    /* calculate jet parameter */
    double t[2*N_tbin],t_s[2*N_tbin],R[2*N_tbin],R_para[2*N_tbin],R_perp[2*N_tbin],gam_j[2*N_tbin],theta_ene_ave[2*N_tbin],n_amb[2*N_tbin];
    double E_j_0=1.0e55*(1.0-cos(theta_j_0));
    double R_j = pow((3.0-k)*E_j_0/2.0/M_PI/A/C/C,1.0/(3.0-k));
    double rfs_max=R_j*r_nd[0],del_ln_rfs=(log(rfs_max)-log(rfs_min))/(double)(N_tbin);
    /* free streaming phase */
    for (i=0;i<N_tbin;i++){
        R[i] = rfs_min*exp(del_ln_rfs*(double)i);
        t[i] = R[i]/C; /* source rest frame */
        t_s[i] = R[i]/gam_j_0/C; /* jet rest frame */
        gam_j[i] = gam_j_0;
        R_perp[i] = R[i]*(2.0*theta_j_0-sin(2.0*theta_j_0))/2.0/M_PI/(1.0-cos(theta_j_0)); /* Eq. (39) */
        R_para[i] = R[i]*sin(theta_j_0)*sin(theta_j_0)/2.0/(1.0-cos(theta_j_0)); /* Eq. (39) */
        theta_ene_ave[i] = (sin(theta_j_0)-theta_j_0*cos(theta_j_0))/(1.0-cos(theta_j_0)); /* Eq. (40) */
        n_amb[i] = A*pow(R_para[i],-k)/M_PRO;
    }
    /* deceleration phase */
    for (i=N_tbin;i<2*N_tbin;i++){
        t[i] = R_j*t_nd[i-N_tbin]/C; /* source rest frame */
        t_s[i] = R_j*t_s_nd[i-N_tbin]/C; /* jet rest frame */
        R[i] = R_j*r_nd[i-N_tbin];
        gam_j[i] = gam_j_dc[i-N_tbin];
        R_perp[i] = R_j*r_nd[i-N_tbin]*(2.0*theta_j_dc[i-N_tbin]-sin(2.0*theta_j_dc[i-N_tbin]))/2.0/M_PI/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (39) */
        R_para[i] = R_j*r_nd[i-N_tbin]*sin(theta_j_dc[i-N_tbin])*sin(theta_j_dc[i-N_tbin])/2.0/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (39) */
        theta_ene_ave[i] = (sin(theta_j_dc[i-N_tbin])-theta_j_dc[i-N_tbin]*cos(theta_j_dc[i-N_tbin]))/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (40) */
        n_amb[i] = A*pow(R_para[i],-k)/M_PRO;
    }
    
    
    /* outputting the microphysics parameters at the shock */
    FILE *op;
    char head_op[256]="conical_jet",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head_op,dat);
    op = fopen(output_file_name,"w+");
    fprintf(op,"# t[s], t_s[s], gam_j, R[cm], R_perp[cm], R_para[cm], theta_j, n_amb[cm^-3]\n");
    for (i=0;i<2*N_tbin;i++){
        fprintf(op,"%le %le %le %le %le %le %le %le \n",
                t[i],t_s[i],gam_j[i],R[i],R_perp[i],R_para[i],theta_ene_ave[i],n_amb[i]);
    }
    fclose(op);
    
}


void calc_shocked_jet()
{
    /* reading the input file of the jet dynamics */
    double gam_j[2*N_tbin],n_amb[2*N_tbin];
    FILE *ip;
    char head_ip[256]="conical_jet",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    int i=0;
    while (fscanf(ip,"%*le %*le %le %*le %*le %*le %*le %le \n",&gam_j[i],&n_amb[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    
    /* calculating the physical parameters at the (forward) shock downstream*/
    double n_f[2*N_tbin],B_f[2*N_tbin],gam_e_inj[2*N_tbin],gam_e_max[2*N_tbin],gam_e_th[2*N_tbin];
    
    for (i=0;i<2*N_tbin;i++){
        n_f[i] = 4.0*gam_j[i]*n_amb[i];
        B_f[i] = sqrt(32.0*M_PI*eps_B*M_PRO*n_amb[i]*gam_j[i]*(gam_j[i]-1.0))*C;
        gam_e_th[i] = XI_T*M_PRO/M_ELE*(gam_j[i]-1.0);
        gam_e_inj[i] = ZETA_E*M_PRO/M_ELE*(gam_j[i]-1.0);
        gam_e_max[i] = sqrt(9.0*M_PI*ELEC/10.0/SIGMA_T/B_f[i]);
    }
    
    
    /* outputting the microphysics parameters at the shock */
    FILE *op;
    char head_op[256]="shock_acc_conical",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head_op,dat);
    op = fopen(output_file_name,"w+");
    fprintf(op,"# n_f[cm^-3], B_f[G], gam_e_th, gam_e_inj, gam_e_max \n");
    for (i=0;i<2*N_tbin;i++){
        fprintf(op,"%le %le %le %le %le \n",
                n_f[i],B_f[i],gam_e_th[i],gam_e_inj[i],gam_e_max[i]);
    }
    fclose(op);
}


void calc_sync_map()
{
    int i=0,j;
    double t[2*N_tbin],t_s[2*N_tbin],gam_j[2*N_tbin],R_para[2*N_tbin],n_f[2*N_tbin],B_f[2*N_tbin],gam_e_th[2*N_tbin],gam_e_inj[2*N_tbin],gam_e_max[2*N_tbin];
    FILE *ip;
    
    char head_ip[256]="conical_jet",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %*le %*le %le %*le %*le \n",
                  &t[i],&t_s[i],&gam_j[i],&R_para[i])!=EOF){
        i++;
    }
    fclose(ip);
    i=0;
    
    char head_ip2[256] ="shock_acc_conical";
    sprintf(input_file_name,"%s%s%s",path,head_ip2,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %le %le \n",
                  &n_f[i],&B_f[i],&gam_e_th[i],&gam_e_inj[i],&gam_e_max[i])!=EOF){
        i++;
    }
    fclose(ip);
    
    /* setting variables */
    double ne[N_ne],gam_e[N_ne],dene_e[N_ne],dne_dt_inj[N_ne],P_cool[N_ne],gam_ph[N_nph],dene_ph[N_nph],P_nu_syn[N_nph],alpha_nu_syn[N_nph];
    double del_ln_gam_e=0.0,del_ln_gam_ph=0.0;
    double dt=0.0,dr=0.0;
    set_integ_base_ne(ne,gam_e,dene_e,&del_ln_gam_e);
    set_integ_base_pnu(gam_ph,dene_ph,&del_ln_gam_ph);
    
    /* calculate emission intensity */
    FILE *op1,*op2;
    char head1[256]="map_ne",head2[256]="map_pnu",output_file_name1[256]={"\0"},output_file_name2[256]={"\0"};
    sprintf(output_file_name1,"%s%s%s",path,head1,dat);
    sprintf(output_file_name2,"%s%s%s",path,head2,dat);
    op1 = fopen(output_file_name1,"w+");
    op2 = fopen(output_file_name2,"w+");
    for (i=0;i<2*N_tbin-1;i++) {
        printf("t = %12.3e [s]\n",t[i]);
        elec_injection(t_s[i],n_f[i],B_f[i],gam_e_inj[i],gam_e_th[i],gam_e_max[i],gam_e,dne_dt_inj);
        elec_cooling(t_s[i],B_f[i],P_cool,gam_e);
        dt = t_s[i+1]-t_s[i];
        dr = (R_para[i+1]-R_para[i])/(4.0*gam_j[i]);
        elec_time_evolution(dt,gam_e,dene_e,ne,ne,dne_dt_inj,P_cool);
        syn_spec(B_f[i],dr,ne,gam_e,dene_e,del_ln_gam_e,gam_ph,P_nu_syn,alpha_nu_syn);
        for (j=0;j<N_ne;j++){
            fprintf(op1,"%le ",ne[j]);
        }
        fprintf(op1,"\n");
        for (j=0; j<N_nph; j++) {
            fprintf(op2,"%le ",P_nu_syn[j]);
        }
        fprintf(op2,"\n");
    }
    fclose(op1);
    fclose(op2);
    
}

void set_integ_base_ne(double ne[], double gam_e[], double dene_e[], double *del_ln_gam_e)
{
    /* electron energy in unit of MeC2 */
    int i;
    double del_ln_gam_e_tmp = (log(GAMMA_ELEC_MAX)-log(GAMMA_ELEC_MIN))/(double)(N_ne-1);
    for (i=0;i<N_ne;i++){
        gam_e[i] = GAMMA_ELEC_MIN*exp(del_ln_gam_e_tmp*(double)i);
        dene_e[i] = gam_e[i]*MeC2*(exp(del_ln_gam_e_tmp)-1.0);
        ne[i] = 0.0;
    }
    *del_ln_gam_e = del_ln_gam_e_tmp;
    
}

void set_integ_base_pnu(double gam_ph[], double dene_ph[], double *del_ln_gam_ph)
{
    /* photon energy in unit of MeC2 */
    int i;
    double del_ln_gam_ph_tmp = (log(ENE_PH_MAX_DIV_MeC2)-log(ENE_PH_MIN_DIV_MeC2))/(double)(N_nph-1);
    for (i=0;i<N_nph;i++){
        gam_ph[i] = ENE_PH_MIN_DIV_MeC2*exp(del_ln_gam_ph_tmp*(double)i);
        dene_ph[i] = gam_ph[i]*MeC2*(exp(del_ln_gam_ph_tmp)-1.0);
    }
    *del_ln_gam_ph = del_ln_gam_ph_tmp;
    
}

void set_integ_base_mu(double mu[], double *del_mu)
{
    /* scattering angle */
    int i;
    double del_mu_tmp = 2.0/(double)(N_mu-1);
    for (i=0;i<N_mu;i++){
        mu[i] = -1.0+del_mu_tmp*(double)i;
    }
    *del_mu = del_mu_tmp;
    
}


void elec_injection(double t_sh, double n_sh, double B_sh, double gam_e_inj, double gam_e_th, double gam_e_max, double gam_e[], double dne_dt_inj[])
{
    double integ_th=0.0,integ_nth=0.0,gam_e_th_tmp=1.0;
    double del_ln_gam_th = log(gam_e_max)/(double)(N_ne-1);
    int i;
    
    for (i=1;i<N_ne-1;i++)
    {
        gam_e_th_tmp = exp(del_ln_gam_th*(double)i);
        integ_th += (gam_e_th_tmp*sqrt(gam_e_th_tmp*gam_e_th_tmp-1.0)*exp(-gam_e_th_tmp/gam_e_th)*gam_e_th_tmp)*del_ln_gam_th;
    }
    integ_th += 0.5*(gam_e_max*sqrt(gam_e_max*gam_e_max-1.0)*exp(-gam_e_max/gam_e_th)*gam_e_max)*del_ln_gam_th;
    
    if (POW_ELE != 1.0){
        integ_nth = (pow(gam_e_inj,1.0-POW_ELE)-pow(gam_e_max,1.0-POW_ELE))/(POW_ELE-1.0);
    } else {
        integ_nth = log(gam_e_max/gam_e_inj);
    }
    
    double dne_tot = n_sh/t_sh;
    double norm_th = (1.0-FRAC_E)*dne_tot/integ_th;
    double norm_nth = FRAC_E*dne_tot/integ_nth;
    
    for (i=0;i<N_ne;i++){
        if (gam_e[i] > gam_e_inj && gam_e[i] < gam_e_max)
            dne_dt_inj[i] = norm_th*gam_e[i]*sqrt(gam_e[i]*gam_e[i]-1.0)*exp(-gam_e[i]/gam_e_th)+norm_nth*pow(gam_e[i],-POW_ELE);
        else
            dne_dt_inj[i] = norm_th*gam_e[i]*sqrt(gam_e[i]*gam_e[i]-1.0)*exp(-gam_e[i]/gam_e_th);
    }
}


double power_ad(double gam_e, double t)
{
    return gam_e*MeC2/t;
}

double power_syn(double gam_e, double B_sh)
{
    // electron synchrotron energy loss rate (see e.g., Eq. 7.13 of Dermer & Menon)
    // double sin2phi = 2.0/3.0; /* averaging pitch angle */
    // double beta_par = 1.0; /* assuming that particles are relativistic */
    
    return 4.0/3.0*C*SIGMA_T*(B_sh*B_sh/8.0/M_PI)*gam_e*gam_e;
}


void elec_cooling(double t, double B_sh, double P_cool[], double gam_e[])
{
    double P_ad=0.0,P_syn=0.0;
    int i;
    for (i=0;i<N_ne;i++) {
        P_ad = power_ad(gam_e[i],t);
        P_syn = power_syn(gam_e[i],B_sh);
        P_cool[i] = P_ad+P_syn;
    }
}

void elec_time_evolution(double dt, double gam_e[], double dene_e[], double ne_old[], double ne_new[], double dne_dt_inj[], double P_cool[])
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


/* syncrotron emission power in the shocked jet rest frame */
void syn_spec(double B, double dr, double ne[], double gam_e[], double dene_e[], double del_ln_gam_e, double gam_ph[], double P_nu_syn[], double alpha_nu_syn[])
{
    int i,j,k;
    double nu,x,sin_alpha,tau_sa;
    double integ=0.0,integ_alpha=0.0;
    
    sin_alpha = 2.0/3.0;
    for (k=0;k<N_nph;k++) {
        integ = 0.0;
        integ_alpha = 0.0;
        nu = gam_ph[k]*MeC2/H;
        for (i=0;i<N_ne;i++) {
            x= (2.0*M_PI*nu)/(3.0*ELEC*gam_e[i]*gam_e[i]*B/2.0/M_ELE/C*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */
            if (i==0 || i==N_ne-1) {
                integ += 0.5*ne[i]*gam_e[i]*del_ln_gam_e*syn_func_fit(x);
                integ_alpha += -0.5*sin_alpha*pow(gam_e[i],2.0)*(-ne[i]/pow(gam_e[i],2.0))/dene_e[i]*syn_func_fit(x)*gam_e[i]*del_ln_gam_e/MeC2;
            } else {
                integ += ne[i]*gam_e[i]*del_ln_gam_e*syn_func_fit(x);
                integ_alpha += -sin_alpha*pow(gam_e[i],2.0)*(ne[i+1]/pow(gam_e[i+1],2.0)-ne[i]/pow(gam_e[i],2.0))/dene_e[i]*syn_func_fit(x)*gam_e[i]*del_ln_gam_e/MeC2;
            }
        }
        P_nu_syn[k] = sqrt(3.0)*pow(ELEC,3.0)*B*sin_alpha/MeC2*integ; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
        alpha_nu_syn[k] = C*C/8.0/M_PI/nu/nu*sqrt(3.0)*pow(ELEC,3.0)*B*integ_alpha; /* Eq. (6.52) of Rybicki & Lightman */
        
        tau_sa = alpha_nu_syn[k]*dr;
        
        if (tau_sa > 1.0e-6)
            P_nu_syn[k] = (1.0-exp(-tau_sa))*P_nu_syn[k]/tau_sa;
        
        integ = 0.0;
        integ_alpha = 0.0;
    }
}

void distances(double z, double *d_h, double *d_c, double *d_a, double *d_l)
{
  int i,n_int=100000;

  double d_h_tmp = C/H_0;
  double d_c_tmp=0.0,d_a_tmp=0.0,d_l_tmp=0.0;

  double z_tmp = 0.0;
  double dz = z/(double)(n_int-1);
  double integ = 0.0;

  for (i=0;i<n_int;i++){
    if (i == 0 || i == n_int-1)
      integ += (1.0/3.0)/sqrt(OMEGA_M*pow(1.0+z_tmp,3.0)+OMEGA_LAMBDA)*dz; 
    else if (i % 2 == 0)
      integ += (2.0/3.0)/sqrt(OMEGA_M*pow(1.0+z_tmp,3.0)+OMEGA_LAMBDA)*dz;
    else 
      integ += (4.0/3.0)/sqrt(OMEGA_M*pow(1.0+z_tmp,3.0)+OMEGA_LAMBDA)*dz;
    z_tmp += dz;
  }

  d_c_tmp = d_h_tmp*integ;
  d_a_tmp = d_c_tmp/(1.0+z);
  d_l_tmp = d_c_tmp*(1.0+z);

  *d_h = d_h_tmp*PC;
  *d_c = d_c_tmp*PC;
  *d_a = d_a_tmp*PC;
  *d_l = d_l_tmp*PC;
}
