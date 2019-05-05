#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physcon.h"
#include "input.h"
#include "funclist.h"


int main()
{
    calc_lightcurve(nuobs,thetaobs);
    
    return 0;
}


void calc_lightcurve(double nuobs, double thetaobs)
{
    /* calculate luminosity distance */
    double d_h,d_c,d_a,d_l;
    distances(z,&d_h,&d_c,&d_a,&d_l);
    
    /* load Pnu */
    int i,j;
    int rows=2*N_tbin-1;
    int columns=N_nph;
    double **Pnu_list = (double**)malloc(rows*sizeof(double*));
    Pnu_list[0]=(double*)malloc(rows*columns*sizeof(double));
    for(i=0;i<rows;i++){
        Pnu_list[i]=Pnu_list[0]+i*columns;
    }
    FILE *ip;
    char head[256]="map_pnu",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head,dat);
    ip = fopen(input_file_name,"r");
    for (i=0; i<2*N_tbin-1; i++) {
        for (j=0; j<N_nph; j++) {
            fscanf(ip,"%le ",&Pnu_list[i][j]);
        }
        fscanf(ip,"\n");
    }
    fclose(ip);
    
    /* load egg shape */
    double del_ln_tobs=(log(tobs_max)-log(tobs_min))/(double)(N_tobsbin-1),tobs,Fnuobs =0.0;
    double mu_integ[2*N_tbin],dmu_integ[2*N_tbin],beam_fac[2*N_tbin],vol_fac[2*N_tbin];
    int time_index_min,time_index_max;
    int nu_index[2*N_tbin];
    double integ=0.;

    /* integrate Pnu along with the egg shape and output the lightcurve */
    FILE *op;
    char head_op[256]="lc",output_file_name[256]={"\0"};
    sprintf(output_file_name,"%s%s_%1.1e%s",path,head_op,nuobs,dat);
    op = fopen(output_file_name,"w+");
    for (i=0; i<N_tobsbin; i++) {
        tobs=tobs_min*exp(del_ln_tobs*(double)i);
        egg_shape(nuobs,tobs,thetaobs,mu_integ,dmu_integ,beam_fac,vol_fac,&time_index_min,&time_index_max,nu_index);
        for (j=time_index_min;j<time_index_max;j++) {
            integ += Pnu_list[j][nu_index[j]]/pow(beam_fac[j],2.0)*vol_fac[j]*dmu_integ[j];
        }
        Fnuobs = (1.+z)/(2.*d_l*d_l)*integ;
        fprintf(op,"%le %le \n",tobs,Fnuobs);
        integ=0.;
    }
    fclose(op);
    free(Pnu_list);
}

void egg_shape(double nuobs, double tobs, double thetaobs, double mu_integ[], double dmu_integ[], double beam_fac[], double vol_fac[], int *time_index_min, int *time_index_max, int nu_index[])
{
    /* load jet dynamics */
    double t[2*N_tbin],gam_j[2*N_tbin],R[2*N_tbin],theta_j[2*N_tbin];
    FILE *ip;
    char head_ip[256]="conical_jet",dat[256]=".dat",input_file_name[256]={"\0"};
    sprintf(input_file_name,"%s%s%s",path,head_ip,dat);
    ip = fopen(input_file_name,"r");
    fscanf(ip,"%*[^\n]");
    int i=0,j;
    while (fscanf(ip,"%le %*le %le %le %*le %*le %le %*le \n",&t[i],&gam_j[i],&R[i],&theta_j[i])!=EOF) {
        i++;
    }
    fclose(ip);
    
    
    /* find the minimum and maximum time indices */
    i = 0;
    int time_index_min_tmp,time_index_max_tmp;
    while ((t[i]+R[i]/C) < tobs/(1.+z)){
        i++;
    }
    time_index_min_tmp = i;
    while ((t[i]-R[i]/C) < tobs/(1.+z)){
        i++;
    }
    time_index_max_tmp = i;
    int N_mu_integ = time_index_max_tmp-time_index_min_tmp;
    for (i=0; i<2*N_tbin; i++) {
        nu_index[i] = 0;
        dmu_integ[i] = 0.;
        if (i >= time_index_min_tmp && i < time_index_max_tmp){
            mu_integ[i] = (t[i]-tobs/(1.+z))*C/R[i];
            beam_fac[i] = gam_j[i]*(1.-sqrt(1.-1./gam_j[i]/gam_j[i])*mu_integ[i]);
            vol_fac[i] = R[i]*R[i]*(R[i+1]-R[i]);
        } else {
            mu_integ[i] = 2.; /* returning nan with acos() */
            beam_fac[i] = 1.;
            vol_fac[i] = 0.;
        }
    }
    *time_index_min = time_index_min_tmp;
    *time_index_max = time_index_max_tmp;
    
    double gam_ph[N_nph],dene_ph[N_nph];
    double del_ln_gam_ph=0.0,dphi=0.0;
    set_integ_base_pnu(gam_ph,dene_ph,&del_ln_gam_ph);
    for (i=time_index_min_tmp; i<time_index_max_tmp; i++) {
        /* find the nu index */
        j=0;
        while (MeC2*gam_ph[j]/H < nuobs*beam_fac[i]) {
            j++;
        }
        nu_index[i] = j;
        
        /* define the edge of the "egg" */
        if (mu_integ[i] < 0.){  /* counter jet */
            if ( M_PI - thetaobs + theta_j[i] < acos(mu_integ[i]) || M_PI - thetaobs - theta_j[i] > acos(mu_integ[i]) ) {
                dphi = 0.;
            } else if (theta_j[i] - thetaobs + acos(mu_integ[i]) < M_PI){
                dphi = 2.*acos((cos(theta_j[i]) + mu_integ[i]*cos(thetaobs))/fabs(sqrt(1.-mu_integ[i]*mu_integ[i]))/sin(thetaobs));
            } else {
                dphi = 2.*M_PI;
            }
        } else {
            if ( thetaobs - theta_j[i] > acos(mu_integ[i]) || thetaobs + theta_j[i] < acos(mu_integ[i]) ) {
                dphi = 0.;
            } else if ( thetaobs + acos(mu_integ[i]) - theta_j[i] > 0. ){
                dphi = 2.*acos((cos(theta_j[i])-mu_integ[i]*cos(thetaobs))/sqrt(1.-mu_integ[i]*mu_integ[i])/sin(thetaobs));
            } else {
                dphi = 2.*M_PI;
            }
        }
        
        if (i == time_index_min_tmp){
            dmu_integ[i] = (mu_integ[i]+1.)*(dphi/2./M_PI);
        } else {
            dmu_integ[i] = (mu_integ[i]-mu_integ[i-1])*(dphi/2./M_PI);
        }
    }
    
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
