#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physcon.h"
#include "input.h"
#include "funclist.h"


int main()
{
    calc_shocked_jet();
    calc_sync_map();
    
    return 0;
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
        gam_e_max[i] = GAMMA_ELEC_MAX;//sqrt(9.0*M_PI*ELEC/10.0/SIGMA_T/B_f[i]);
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
    double t[2*N_tbin],t_j[2*N_tbin],gam_j[2*N_tbin],R_para[2*N_tbin],n_f[2*N_tbin],B_f[2*N_tbin],gam_e_th[2*N_tbin],gam_e_inj[2*N_tbin],gam_e_max[2*N_tbin];
    FILE *ip;
    
    char head_ip[256]="conical_jet",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    while (fscanf(ip,"%le %le %le %*le %*le %le %*le %*le \n",
                  &t[i],&t_j[i],&gam_j[i],&R_para[i])!=EOF){
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
        elec_injection(t_j[i],n_f[i],B_f[i],gam_e_inj[i],gam_e_th[i],gam_e_max[i],gam_e,dne_dt_inj);
        elec_cooling(t_j[i],B_f[i],P_cool,gam_e);
        dt = t_j[i+1]-t_j[i];
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
        
        if (tau_sa > 1.0e-6){
            P_nu_syn[k] = (1.0-exp(-tau_sa))*P_nu_syn[k]/tau_sa;
        }
        
        integ = 0.0;
        integ_alpha = 0.0;
    }
}
