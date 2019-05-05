#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physcon.h"
#include "input.h"
#include "funclist.h"


int main()
{
    calc_conical_model();
    calc_jet();
    
    return 0;
}


void calc_conical_model()
{
    FILE *op;
    char head[256]="conical_model",dat[256]=".dat",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head,dat);
    
    int i;
    double r_s_min=pow((1.0-cos(theta_j_0))*(gam_j_0*gam_j_0-1.0),-1.0/(3.0-k)),r_s_max=1.0e1,del_ln_r_s=(log(r_s_max)-log(r_s_min))/(double)(N_tbin-1),r_s=r_s_min,dr_s=0.0;
    double theta_j_tmp=theta_j_0,u_tmp=sqrt(gam_j_0*gam_j_0-1.0);
    double t=r_s_min,dt=0.,t_j=t/(sqrt(u_tmp*u_tmp+1.));
    
    op = fopen(output_file_name,"w+");
    fprintf(op,"#k=%le A=%le a=%le b=%le\n",k,A,a,b);
    for (i=0;i<N_tbin;i++){
        r_s = r_s_min*exp(del_ln_r_s*(double)i);
        dr_s = r_s*(exp(del_ln_r_s)-1.0);
        dt = dr_s*sqrt(u_tmp*u_tmp+1.)/u_tmp;
        t += dt; /* source rest frame */
        t_j += dt/(sqrt(u_tmp*u_tmp+1.)); /* jet rest frame */
        
        theta_j_tmp += 1.0/r_s/pow(1.0+pow(r_s,k-3.0)/(1.0-cos(theta_j_tmp)),(1.0+a)/2.0)/pow(theta_j_tmp,a)*dr_s;
        if (theta_j_tmp>M_PI/2.0)
            theta_j_tmp = M_PI/2.0;
        u_tmp = pow(r_s,-(3.0-k)/2.0)/sqrt(1.0-cos(theta_j_tmp));
        
        fprintf(op,"%le %le %le %le %le %le \n",
                r_s*pow(4.0/(3.0-k),-1.0/(3.0-k)),t*pow(4.0/(3.0-k),-1.0/(3.0-k)),t_j*pow(4.0/(3.0-k),-1.0/(3.0-k)),theta_j_tmp,u_tmp,sqrt(u_tmp*u_tmp+1.0));
    }
    fclose(op);
}


void calc_jet()
{
    /* reading the input file of the jet dynamics */
    double r_nd[N_tbin],t_nd[N_tbin],t_j_nd[N_tbin],theta_j_dc[N_tbin],gam_j_dc[N_tbin];
    FILE *ip;
    char head_ip[256]="conical_model",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    int i=0;
    while (fscanf(ip,"%le %le %le %le %*le %le \n",&r_nd[i],&t_nd[i],&t_j_nd[i],&theta_j_dc[i],&gam_j_dc[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    
    /* calculate jet parameter */
    double t[2*N_tbin],t_j[2*N_tbin],R[2*N_tbin],R_para[2*N_tbin],R_perp[2*N_tbin],gam_j[2*N_tbin],theta_ene_ave[2*N_tbin],n_amb[2*N_tbin];
    double E_j_0=1.0e55*(1.0-cos(theta_j_0));
    double R_j = pow((3.0-k)*E_j_0/2.0/M_PI/A/C/C,1.0/(3.0-k));
    double rfs_max=R_j*r_nd[0],del_ln_rfs=(log(rfs_max)-log(rfs_min))/(double)(N_tbin);
    /* free streaming phase */
    for (i=0;i<N_tbin;i++){
        R[i] = rfs_min*exp(del_ln_rfs*(double)i);
        t[i] = R[i]/C; /* source rest frame */
        t_j[i] = R[i]/gam_j_0/C; /* jet rest frame */
        gam_j[i] = gam_j_0;
        R_perp[i] = R[i]*(2.0*theta_j_0-sin(2.0*theta_j_0))/2.0/M_PI/(1.0-cos(theta_j_0)); /* Eq. (39) */
        R_para[i] = R[i]*sin(theta_j_0)*sin(theta_j_0)/2.0/(1.0-cos(theta_j_0)); /* Eq. (39) */
        theta_ene_ave[i] = (sin(theta_j_0)-theta_j_0*cos(theta_j_0))/(1.0-cos(theta_j_0));
        n_amb[i] = A*pow(R[i],-k)/M_PRO;
    }
    /* deceleration phase */
    for (i=N_tbin;i<2*N_tbin;i++){
        t[i] = R_j*t_nd[i-N_tbin]/C; /* source rest frame */
        t_j[i] = R_j*t_j_nd[i-N_tbin]/C; /* jet rest frame */
        R[i] = R_j*r_nd[i-N_tbin];
        gam_j[i] = gam_j_dc[i-N_tbin];
        R_perp[i] = R_j*r_nd[i-N_tbin]*(2.0*theta_j_dc[i-N_tbin]-sin(2.0*theta_j_dc[i-N_tbin]))/2.0/M_PI/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (39) */
        R_para[i] = R_j*r_nd[i-N_tbin]*sin(theta_j_dc[i-N_tbin])*sin(theta_j_dc[i-N_tbin])/2.0/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (39) */
        theta_ene_ave[i] = (sin(theta_j_dc[i-N_tbin])-theta_j_dc[i-N_tbin]*cos(theta_j_dc[i-N_tbin]))/(1.0-cos(theta_j_dc[i-N_tbin])); /* Eq. (40) */
        n_amb[i] = A*pow(R[i],-k)/M_PRO;
    }
    
    
    /* outputting the microphysics parameters at the shock */
    FILE *op;
    char head_op[256]="conical_jet",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s%s",path,head_op,dat);
    op = fopen(output_file_name,"w+");
    fprintf(op,"# t[s], t_j[s], gam_j, R[cm], R_perp[cm], R_para[cm], theta_j, n_amb[cm^-3]\n");
    for (i=0;i<2*N_tbin;i++){
        fprintf(op,"%le %le %le %le %le %le %le %le \n",
                t[i],t_j[i],gam_j[i],R[i],R_perp[i],R_para[i],theta_ene_ave[i],n_amb[i]);
    }
    fclose(op);
    
}
