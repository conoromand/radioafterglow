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

const int N_tbin = 1000;
const int N_ne = 512;
const int N_nph = 512; /* need to be the same as N_ne */
const int N_mu = 32;

const double GAMMA_ELEC_MIN = 1.0;/* minimum Lorentz factor of electron */
const double GAMMA_ELEC_MAX = 2.0e9;/* maximum Lorentz factor of electron */
const double ENE_PH_MIN_DIV_MeC2 = 1.0e-14;/* mimimum energy of photon in unit of electron mass */
const double ENE_PH_MAX_DIV_MeC2 = 1.0e3;/* maximum energy of photon in unit of electron mass */

const double ZETA_E = 0.1;//0.4;
const double FRAC_E = 1.0;//0.01;
const double XI_T = 0.24;
const double POW_ELE = 2.0;


void calc_conical_model(char path[256], double k, double A, double a, double b, double theta_j_0, double gam_j_0);
void calc_shocked_jet(char path[256], double k, double A, double theta_j_0, double E_j_0, double eps_B);
void set_var(double mu[], double *del_mu, double ne[], double gam_e[], double dene_e[], double *del_ln_gam_e, double gam_ph[], double dene_ph[], double *del_ln_gam_ph);
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
    /* input parameters */
    /////////////////////////////////////////
    /* Ambient density profile */
    double v_w_str=2.0e8;
    double Mdot_str=1.0e-4*M_SUN/YR;
    double k=2.0;
    double A=Mdot_str/4.0/M_PI/v_w_str;
    //double k=0.0;
    //double A=0.1*M_PRO;
    double a=1.0; /* see Granot & Piran 2011 */
    double b=0.45; /* see Granot & Piran 2011 */
    /* redshift */
    double z = 0.1;
    /* jet parameters  */
    double theta_j_0=2.0*M_PI*5.0/360.0;
    double gam_j_0=200.0;
    double E_j_0=1.0e55*(1.0-cos(theta_j_0));
    /* shock parameters */
    double eps_B = 0.001;
    double eps_e = 0.01;
    double p = 2.2;
    /* path */
    char path[256]="/Users/kakashi/offaxis_data/";
    /////////////////////////////////////////
    
    
    /* setting variables */
    /////////////////////////////////////////
    double mu[N_mu],ne[N_ne],gam_e[N_ne],dene_e[N_ne],dne_dt_inj[N_ne],P_cool[N_ne],gam_ph[N_nph],dene_ph[N_nph],P_nu_syn[N_nph],alpha_nu_syn[N_nph];
    double del_mu=0.0,del_ln_gam_e=0.0,del_ln_gam_ph=0.0;
    double dt=0.0,dr=0.0;
    set_var(mu,&del_mu,ne,gam_e,dene_e,&del_ln_gam_e,gam_ph,dene_ph,&del_ln_gam_ph);
    /////////////////////////////////////////

    
    calc_conical_model(path,k,A,a,b,theta_j_0,gam_j_0);
    calc_shocked_jet(path,k,A,theta_j_0,E_j_0,eps_B);
    
    
    int i=0,j;
    double t_s[N_tbin],t_prime[N_tbin],gam_j[N_tbin],R_perp[N_tbin],R_para[N_tbin],theta_ene_ave[N_tbin],n_f[N_tbin],B_f[N_tbin],gam_e_th[N_tbin],gam_e_inj[N_tbin],gam_e_max[N_tbin];
    FILE *ip;
    char head_ip[256]="shock_acc_conical",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %le %le %le %le %le %le %le %le \n",
           &t_s[i],&t_prime[i],&gam_j[i],&R_perp[i],&R_para[i],&theta_ene_ave[i],&n_f[i],&B_f[i],&gam_e_th[i],&gam_e_inj[i],&gam_e_max[i])!=EOF) {
        i++;
    }
    fclose(ip);

    /* calculate emission intensity */
    FILE *op1,*op2;
    char head1[256]="map_ne",head2[256]="map_pnu",output_file_name1[256]={"\0"},output_file_name2[256]={"\0"};
    sprintf(output_file_name1,"%s%s%s",path,head1,dat);
    sprintf(output_file_name2,"%s%s%s",path,head2,dat);
    op1 = fopen(output_file_name1,"w+");
    op2 = fopen(output_file_name2,"w+");
    for (i=0;i<N_tbin-1;i++) {
        printf("t = %12.3e [s]\n",t_s[i]);
        elec_injection(t_prime[i],n_f[i],B_f[i],gam_e_inj[i],gam_e_th[i],gam_e_max[i],gam_e,dne_dt_inj);
        elec_cooling(t_prime[i],B_f[i],P_cool,gam_e);
        dt = t_prime[i+1]-t_prime[i];
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
    
    
    double d_h=0.0,d_c=0.0,d_a=0.0,d_l=0.0;
    distances(z,&d_h,&d_c,&d_a,&d_l);
    
    return 0;
}


void calc_conical_model(char path[256], double k, double A, double a, double b, double theta_j_0, double gam_j_0)
{
    FILE *op;
    char head[256]="conical_model",dat[256]=".dat",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head,dat);
    
    int i;
    double r_s_min=pow((1.0-cos(theta_j_0))*(gam_j_0*gam_j_0-1.0),-1.0/(3.0-k)),r_s_max=1.0e1,del_ln_r_s=(log(r_s_max)-log(r_s_min))/(double)(N_tbin-1),r_s=r_s_min,dr_s=0.0;
    double theta_j_tmp=theta_j_0,u_tmp=sqrt(gam_j_0*gam_j_0-1.0);
    double t=r_s_min,dt=0.0,t_s=(1.0-sqrt(1.0-1.0/(u_tmp*u_tmp+1.0)))*t;
    
    op = fopen(output_file_name,"w+");
    fprintf(op,"#k=%le A=%le a=%le b=%le\n",k,A,a,b);
    for (i=0;i<N_tbin;i++){
        r_s = r_s_min*exp(del_ln_r_s*(double)i);
        dr_s = r_s*(exp(del_ln_r_s)-1.0);
        dt = dr_s/sqrt(1.0-1.0/(u_tmp*u_tmp+1.0));
        t += dt;
        t_s += (1.0-sqrt(1.0-1.0/(u_tmp*u_tmp+1.0)))*dt;
        
        theta_j_tmp += 1.0/r_s/pow(1.0+pow(r_s,k-3.0)/(1.0-cos(theta_j_tmp)),(1.0+a)/2.0)/pow(theta_j_tmp,a)*dr_s;
        if (theta_j_tmp>M_PI/2.0)
            theta_j_tmp = M_PI/2.0;
        u_tmp = pow(r_s,-(3.0-k)/2.0)/sqrt(1.0-cos(theta_j_tmp));
        
        fprintf(op,"%le %le %le %le %le %le \n",
                r_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),t*pow(4.0/(3.0-k),-1.0/(3.0-k)),t_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),theta_j_tmp,u_tmp,sqrt(u_tmp*u_tmp+1.0));
    }
    fclose(op);
}


void calc_shocked_jet(char path[256], double k, double A, double theta_j_0, double E_j_0, double eps_B)
{
    int i=0;
    
    /* reading the input file of the jet dynamics */
    double r[N_tbin],t[N_tbin],t_s[N_tbin],t_prime[N_tbin],theta_j[N_tbin],gam_j[N_tbin];
    FILE *ip;
    char head_ip[256]="conical_model",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %le %*le %le \n",&r[i],&t[i],&t_s[i],&theta_j[i],&gam_j[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    double R_para[N_tbin],R_perp[N_tbin],theta_ene_ave[N_tbin],n_amb[N_tbin];
    double R_j = pow((3.0-k)*E_j_0/2.0/M_PI/A/C/C,1.0/(3.0-k));
    
    for (i=0;i<N_tbin;i++){
        t[i] = R_j*t[i]/C;
        t_s[i] = R_j*t_s[i]/C; /* source rest frame */
        t_prime[i] = t_s[i]/gam_j[i]; /* jet rest frame */
        n_amb[i] = A*pow(R_para[i],-k)/M_PRO;
        R_perp[i] = R_j*r[i]*(2.0*theta_j[i]-sin(2.0*theta_j[i]))/2.0/M_PI/(1.0-cos(theta_j[i]));
        R_para[i] = R_j*r[i]*sin(theta_j[i])*sin(theta_j[i])/2.0/(1.0-cos(theta_j[i]));
        theta_ene_ave[i] = (sin(theta_j[i])-theta_j[i]*cos(theta_j[i]))/(1.0-cos(theta_j[i]));
    }
    
    /* calculating the physical parameters at the (forward) shock downstream*/
    double n_f[N_tbin],B_f[N_tbin],gam_e_inj[N_tbin],gam_e_max[N_tbin],gam_e_th[N_tbin];
    
    for (i=0;i<N_tbin;i++){
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
    fprintf(op,"# t_s[s], t_prime[s], gam_j, R_perp[cm], R_para[cm], theta_j, n_f[cm^-3], B_f[G], gam_e_th, gam_e_inj, gam_e_max \n");
    for (i=0;i<N_tbin;i++){
        fprintf(op,"%le %le %le %le %le %le %le %le %le %le %le \n",
                t_s[i],t_prime[i],gam_j[i],R_perp[i],R_para[i],theta_ene_ave[i],n_f[i],B_f[i],gam_e_th[i],gam_e_inj[i],gam_e_max[i]);
    }
    fclose(op);
}


void set_var(double mu[], double *del_mu, double ne[], double gam_e[], double dene_e[], double *del_ln_gam_e, double gam_ph[], double dene_ph[], double *del_ln_gam_ph)
{
    /* setting basic variables as strings */
    
    /* electron energy in unit of MeC2 */
    int i;
    double del_ln_gam_e_tmp = (log(GAMMA_ELEC_MAX)-log(GAMMA_ELEC_MIN))/(double)(N_ne-1);
    for (i=0;i<N_ne;i++){
        gam_e[i] = GAMMA_ELEC_MIN*exp(del_ln_gam_e_tmp*(double)i);
        dene_e[i] = gam_e[i]*MeC2*(exp(del_ln_gam_e_tmp)-1.0);
        ne[i] = 0.0;
    }
    *del_ln_gam_e = del_ln_gam_e_tmp;
    
    /* photon energy in unit of MeC2 */
    double del_ln_gam_ph_tmp = (log(ENE_PH_MAX_DIV_MeC2)-log(ENE_PH_MIN_DIV_MeC2))/(double)(N_nph-1);
    for (i=0;i<N_nph;i++){
        gam_ph[i] = ENE_PH_MIN_DIV_MeC2*exp(del_ln_gam_ph_tmp*(double)i);
        dene_ph[i] = gam_ph[i]*MeC2*(exp(del_ln_gam_ph_tmp)-1.0);
    }
    *del_ln_gam_ph = del_ln_gam_ph_tmp;
    
    /* scattering angle */
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
