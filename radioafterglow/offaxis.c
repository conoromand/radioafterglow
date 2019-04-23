#include <stdio.h>
#include <math.h>

const int n_rbin = 10000;

const double E_c = 4.3e54; /* [erg] */
const double M_c = 22.0*2e33; /* [g] */
const double R_c = 2.0e13; /* [cm] */
const double kappa = 0.35; /* [cm^2 g^{-1}] */
const double T_ion = 6000; /* [K] */
const double z = 3.0; /* redshift */

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

const int N_ne = 256;
const int N_nph = 256; /* need to be the same as N_ne */

const double XI_T = 0.24;
const double ZETA_E = 0.4;


double theta_ana(double r, double a, double b, double k);
double gam_ana(double r, double a, double b, double k);
void conical_model(double k, double A, double a, double b, double z, double theta_j_0, double gam_j_0);
void shocked_jet_conical(double k, double A, double z, double theta_j_0, double E_j_0, double eps_B);
void sync_spec_para_conical(double z, double d_l, double k, double p);
double sync_spec_conical(double nu, double nu_m, double nu_c, double F_nu_max, double gam_j, double B_f, double R_perp, double p, double z, double d_l);
double syn_func_fit(double x);
void syn_spec(double B, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double mu[], double del_mu, double gamma_ph[], double P_nu_syn[], double alpha_nu_syn[]);
void distances(double z, double *d_h, double *d_c, double *d_a, double *d_l);

int main()
{
  /* Ambient density profile */
  double v_w_str=2.0e8;
  double Mdot_str=1.0e-4*M_SUN/YR;
  double k=2.0;
  double A=Mdot_str/4.0/M_PI/v_w_str;
  //double k=0.0;
  //double A=0.1*M_PRO;

  /* Granot & Piran 2011 */
  double a=1.0; 
  double b=0.45; 

  /* redshift and distance */
  double z = 20.0;
  double d_h=0.0,d_c=0.0,d_a=0.0,d_l=0.0;
  distances(z,&d_h,&d_c,&d_a,&d_l);

  /* jet parameters  */
  double theta_j_0=2.0*M_PI*5.0/360.0;
  double gam_j_0=200.0;
  double E_j_0=1.0e55*(1.0-cos(theta_j_0));

  /* shock parameters */
  double eps_B = 0.001;
  double eps_e = 0.01;
  double p = 2.2;

  conical_model(k,A,a,b,z,theta_j_0,gam_j_0);
  shocked_jet_conical(k,A,z,theta_j_0,E_j_0,eps_B);
    printf("test\n");

  return 0;
}

double theta_ana(double r, double a, double b, double k)
{
  /* Eq.(25) of Granot & Piran 2011 */
  return b*pow(r,-(3.0-k)*(1.0+a)/(3.0+a)/(3.0+a))*exp((3.0+a)/(3.0-k)/(1.0+a)*pow(r,(3.0-k)*(1.0+a)/(3.0+a)));
}

double gam_ana(double r, double a, double b, double k)
{
  /* Eq.(26) of Granot & Piran 2011 */
  return 1.0/b*pow(r,-2.0*(3.0-k)/(3.0+a)/(3.0+a))*exp(-(3.0+a)/(3.0-k)/(1.0+a)*pow(r,(3.0-k)*(1.0+a)/(3.0+a)));
}


void conical_model(double k, double A, double a, double b, double z, double theta_j_0, double gam_j_0)
{
  FILE *op;
  char head[256]="conical_model",dat[256]=".dat",output_file_name[256]={"\0"};
  sprintf(output_file_name,"%s%s",head,dat);

  int i;
  double r_s_min=pow((1.0-cos(theta_j_0))*(gam_j_0*gam_j_0-1.0),-1.0/(3.0-k)),r_s_max=1.0e1,del_ln_r_s=(log(r_s_max)-log(r_s_min))/(double)(n_rbin-1),r_s=r_s_min,dr_s=0.0;
  double theta_j_tmp=theta_j_0,u_tmp=sqrt(gam_j_0*gam_j_0-1.0);
  double t=r_s_min,dt=0.0,t_obs=(1.0+z)*(1.0-sqrt(1.0-1.0/(u_tmp*u_tmp+1.0)))*t;

  op = fopen(output_file_name,"w+");
  fprintf(op,"#k=%le A=%le a=%le b=%le\n",k,A,a,b);
  for (i=0;i<n_rbin;i++){
    r_s = r_s_min*exp(del_ln_r_s*(double)i);
    dr_s = r_s*(exp(del_ln_r_s)-1.0);
    dt = dr_s/sqrt(1.0-1.0/(u_tmp*u_tmp+1.0));
    t += dt;
    t_obs += (1.0+z)*(1.0-sqrt(1.0-1.0/(u_tmp*u_tmp+1.0)))*dt;

    theta_j_tmp += 1.0/r_s/pow(1.0+pow(r_s,k-3.0)/(1.0-cos(theta_j_tmp)),(1.0+a)/2.0)/pow(theta_j_tmp,a)*dr_s;
    if (theta_j_tmp>M_PI/2.0)
      theta_j_tmp = M_PI/2.0;
    u_tmp = pow(r_s,-(3.0-k)/2.0)/sqrt(1.0-cos(theta_j_tmp));

    fprintf(op,"%le %le %le %le %le %le \n",
	    r_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),t*pow(4.0/(3.0-k),-1.0/(3.0-k)),t_obs*pow(4.0/(3.0-k),-1.0/(3.0-k)),theta_j_tmp,u_tmp,sqrt(u_tmp*u_tmp+1.0));
  }
  fclose(op);
}


void shocked_jet_conical(double k, double A, double z, double theta_j_0, double E_j_0, double eps_B)
{
    int i=0;
    
    /* reading the input file of the jet dynamics */
    double r[n_rbin],t[n_rbin],t_obs[n_rbin],theta_j[n_rbin],gam_j[n_rbin];
    FILE *ip;
    char head_ip[256]="conical_model",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s",head_ip,dat);
    
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %le %*le %le \n",&r[i],&t[i],&t_obs[i],&theta_j[i],&gam_j[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    /* calculating the microphysics parameters at the shock */
    double R_para[n_rbin],R_perp[n_rbin],theta_ene_ave[n_rbin];
    double n_amb,n_f[n_rbin],B_f[n_rbin],gam_e_inj[n_rbin],gam_e_max[n_rbin],gam_e_th[n_rbin];
    double R_j = pow((3.0-k)*E_j_0/2.0/M_PI/A/C/C,1.0/(3.0-k));
    
    for (i=0;i<n_rbin;i++){
        t[i] = R_j*t[i]/C;
        t_obs[i] = R_j*t_obs[i]/C;
        R_perp[i] = R_j*r[i]*(2.0*theta_j[i]-sin(2.0*theta_j[i]))/2.0/M_PI/(1.0-cos(theta_j[i]));
        R_para[i] = R_j*r[i]*sin(theta_j[i])*sin(theta_j[i])/2.0/(1.0-cos(theta_j[i]));
        theta_ene_ave[i] = (sin(theta_j[i])-theta_j[i]*cos(theta_j[i]))/(1.0-cos(theta_j[i]));
        n_amb = A*pow(R_para[i],-k)/M_PRO;
        n_f[i] = 4.0*gam_j[i]*n_amb;
        B_f[i] = sqrt(32.0*M_PI*eps_B*M_PRO*n_amb*gam_j[i]*(gam_j[i]-1.0))*C;
        gam_e_th[i] = XI_T*M_PRO/M_ELE*(gam_j[i]-1.0);
        gam_e_inj[i] = ZETA_E*M_PRO/M_ELE*(gam_j[i]-1.0);
        gam_e_max[i] = sqrt(9.0*M_PI*ELEC/10.0/SIGMA_T/B_f[i]);
    }
    
    /* outputting the microphysics parameters at the shock */
    FILE *op;
    char head_op[256]="shock_acc_conical",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s",head_op,dat);
    op = fopen(output_file_name,"w+");
    fprintf(op,"# t[s], t_obs[s], gam_j, R_perp[cm], R_para[cm], theta_j, n_f[cm^-3], B_f[G], gam_e_th, gam_e_inj, gam_e_max \n");
    for (i=0;i<n_rbin;i++){
        fprintf(op,"%le %le %le %le %le %le %le %le %le %le %le \n",
                t[i],t_obs[i],gam_j[i],R_perp[i],R_para[i],theta_ene_ave[i],n_f[i],B_f[i],gam_e_th[i],gam_e_inj[i],gam_e_max[i]);
    }
    fclose(op);
    
}

/* energy distribution of electrons in the shocked jet rest frame */
void elec_spec(double n_sh){
    
}

/* integration of Eqs. (2) and (4) of Granot, Piran, and Sari 99 */
/* 先ずは愚直に1000time step分のP_nuとalpha_nuを吐いてみよう。*/

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
void syn_spec(double B, double ne[], double gamma_e[], double dene_e[], double del_ln_gamma_e, double mu[], double del_mu, double gamma_ph[], double P_nu_syn[], double alpha_nu_syn[])
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
            x= (2.0*M_PI*nu)/(3.0*ELEC*gamma_e[i]*gamma_e[i]*B/2.0/M_ELE/C*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */
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
        P_nu_syn[k] = sqrt(3.0)*pow(ELEC,3.0)*B/MeC2*integ; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
        alpha_nu_syn[k] = C*C/8.0/M_PI/nu/nu*sqrt(3.0)*pow(ELEC,3.0)*B*integ_alpha; /* Eq. (6.52) of Rybicki & Lightman */
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














void sync_spec_para_conical(double z, double d_l, double k, double p)
{
    int i=0;
    
    /* reading the input file of the shock accreleation */
    double t[n_rbin],t_obs[n_rbin],gam_j[n_rbin],R_perp[n_rbin],R_para[n_rbin],theta_ene_ave[n_rbin],n_amb[n_rbin],B_f[n_rbin],gam_e_m[n_rbin],gam_e_c[n_rbin];
    
    FILE *ip;
    char head_ip[256]="shock_acc_conical",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s",head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %le %le %le %le %le %le %le \n",
                  &t[i],&t_obs[i],&gam_j[i],&R_perp[i],&R_para[i],&theta_ene_ave[i],&n_amb[i],&B_f[i],&gam_e_m[i],&gam_e_c[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    /* calculating the synchrotron spectrum */
    double F_nu_max[n_rbin],nu_m[n_rbin],nu_c[n_rbin];
    for (i=0;i<n_rbin;i++){
        F_nu_max[i] = SIGMA_T*M_ELE*C*C*pow(R_para[i],3.0)*n_amb[i]*B_f[i]*gam_j[i]*(1.0+z)/3.0/(3.0-k)/ELEC/d_l/d_l;
        nu_m[i] = gam_j[i]*gam_e_m[i]*gam_e_m[i]*ELEC*B_f[i]/2.0/M_PI/M_ELE/C/(1.0+z);
        nu_c[i] = gam_j[i]*gam_e_c[i]*gam_e_c[i]*ELEC*B_f[i]/2.0/M_PI/M_ELE/C/(1.0+z);
    }
    
    /* outputting the microphysics parameters at the shock */
    FILE *op;
    char head_op[256]="sync_lc_conical",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s",head_op,dat);
    op = fopen(output_file_name,"w+");
    fprintf(op,"# t_obs[s], nu_m[Hz], nu_c[c], F_nu_max[erg/s/cm^2/Hz] \n");
    for (i=0;i<n_rbin;i++){
        fprintf(op,"%le %le %le %le \n",t_obs[i],nu_m[i],nu_c[i],F_nu_max[i]);
    }
    fclose(op);
}

double sync_spec_conical(double nu, double nu_m, double nu_c, double F_nu_max, double gam_j, double B_f, double R_perp, double p, double z, double d_l)
{
    double gam_e_char_obs,F_nu_BB_crit,F_nu;
    gam_e_char_obs = sqrt(2.0*M_PI*M_ELE*C*(1.0+z)*nu/gam_j/ELEC/B_f);
    F_nu_BB_crit = 2.0*M_PI*nu*nu*gam_j*gam_e_char_obs*M_ELE*pow(R_perp/d_l,2.0)*pow(1.0+z,3.0);
    
    if (nu_m <= nu_c){
        
        if (nu <= nu_m)
        F_nu = pow(nu/nu_m,1.0/3.0)*F_nu_max;
        else if (nu <= nu_c)
        F_nu = pow(nu/nu_m,-(p-1.0)/2.0)*F_nu_max;
        else
        F_nu = pow(nu_c/nu_m,-(p-1.0)/2.0)*pow(nu/nu_c,-p/2.0)*F_nu_max;
        
        if (F_nu >= F_nu_BB_crit)
        F_nu = F_nu_BB_crit;
        
    } else {
        
        if (nu <= nu_c)
        F_nu = pow(nu/nu_c,1.0/3.0)*F_nu_max;
        else if (nu <= nu_m)
        F_nu = pow(nu/nu_c,-1.0/2.0)*F_nu_max;
        else
        F_nu = pow(nu_m/nu_c,-1.0/2.0)*pow(nu/nu_m,-p/2.0)*F_nu_max;
        
        if (F_nu >= F_nu_BB_crit)
        F_nu = F_nu_BB_crit;
        
    }
    
    return F_nu;
}
