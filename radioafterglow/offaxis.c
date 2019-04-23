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

const double H_0 = 67.74*1.0e5/1.0e6; /* [cm/s/pc] */ 
const double OMEGA_M = 0.3089;
const double OMEGA_LAMBDA = 0.6911;
const double OMEGA_IGM = 0.04374;

const double M_SUN = 1.989e33; /* [g]  */
const double PC = 3.086e18; /* [cm] */
const double YR = 365.0*24.0*60.0*60.0; /* [s]  */

const double lam_obs_U = 0.3652e-4; /* [cm] */
const double lam_obs_V = 0.4448e-4; /* [cm] */
const double lam_obs_B = 0.5505e-4; /* [cm] */
const double lam_obs_R = 0.6588e-4; /* [cm] */
const double lam_obs_H = 1.654e-4; /* [cm] */
const double lam_obs_K = 2.179e-4; /* [cm] */
const double lam_obs_L = 3.547e-4; /* [cm] */


double theta_ana(double r, double a, double b, double k);
double gam_ana(double r, double a, double b, double k);
void conical_model(double k, double A, double a, double b, double z, double theta_j_0, double gam_j_0);
void shock_acc_conical(double k, double A, double z, double theta_j_0, double E_j_0, double eps_B, double eps_e, double p);
void sync_specpara_conical(double z, double d_l, double k, double p);
double sync_spec_conical(double nu, double nu_m, double nu_c, double F_nu_max, double gam_j, double B_f, double R_perp, double p, double z, double d_l);
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
  shock_acc_conical(k,A,z,theta_j_0,E_j_0,eps_B,eps_e,p);
  sync_specpara_conical(z,d_l,k,p);

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


void shock_acc_conical(double k, double A, double z, double theta_j_0, double E_j_0, double eps_B, double eps_e, double p)
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
  double n_amb[n_rbin],B_f[n_rbin],gam_e_m[n_rbin],gam_e_c[n_rbin];
  double R_j = pow((3.0-k)*E_j_0/2.0/M_PI/A/C/C,1.0/(3.0-k));

  for (i=0;i<n_rbin;i++){
    t[i] = R_j*t[i]/C;
    t_obs[i] = R_j*t_obs[i]/C;
    R_perp[i] = R_j*r[i]*(2.0*theta_j[i]-sin(2.0*theta_j[i]))/2.0/M_PI/(1.0-cos(theta_j[i]));
    R_para[i] = R_j*r[i]*sin(theta_j[i])*sin(theta_j[i])/2.0/(1.0-cos(theta_j[i]));
    theta_ene_ave[i] = (sin(theta_j[i])-theta_j[i]*cos(theta_j[i]))/(1.0-cos(theta_j[i]));
    n_amb[i] = A*pow(R_para[i],-k)/M_PRO;
    B_f[i] = sqrt(32.0*M_PI*eps_B*M_PRO*n_amb[i]*gam_j[i]*(gam_j[i]-1.0))*C;
    gam_e_m[i] = eps_e*(p-2.0)/(p-1.0)*M_PRO/M_ELE*(gam_j[i]-1.0);
    gam_e_c[i] = 6.0*M_PI*M_ELE*C/SIGMA_T/gam_j[i]/B_f[i]/B_f[i]/t_obs[i]*(1.0+z);
  }

  /* outputting the microphysics parameters at the shock */
  FILE *op;
  char head_op[256]="shock_acc_conical",output_file_name[256]={"\0"};
  sprintf(output_file_name,"%s%s",head_op,dat);
  op = fopen(output_file_name,"w+");
  fprintf(op,"# t[s], t_obs[s], gam_j, R_perp[cm], R_para[cm], theta_j, n_amb[cm^-3], B_f[G], gam_e_m, gam_e_c \n");
  for (i=0;i<n_rbin;i++){
    fprintf(op,"%le %le %le %le %le %le %le %le %le %le \n",
	    t[i],t_obs[i],gam_j[i],R_perp[i],R_para[i],theta_ene_ave[i],n_amb[i],B_f[i],gam_e_m[i],gam_e_c[i]);
  }
  fclose(op);
}

void sync_specpara_conical(double z, double d_l, double k, double p)
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
