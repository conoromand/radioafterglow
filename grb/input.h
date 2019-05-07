#ifndef INPUT_H
#define INPUT_H


/////////////////////////////////////////
/* input parameter */
/////////////////////////////////////////
/* path */
const char path[256]="/Users/kakashi/offaxis_data/test_offaxis/";
//const char path[256]="/Users/kakashi/offaxis_data/cnst0/G200thetaj5/";
//const char path[256]="/Users/kakashi/offaxis_data/wind-5/G200thetaj5/";

/* frequency in the observer frame */
const double nuobs = 1.0e8;
/* observer angle */
const double thetaobs = 2.0*M_PI*15./360.0;
/* redshift */
const double z = 0.1;

/* Ambient density profile */
const double rfs_min = 1.0e12; /* minimum radius of the free expansion phase in unit of [cm] */
const double v_w_str=1.0e8;
const double Mdot_str=1.0e-5*M_SUN/YR;
//const double k=2.0;
//const double A=Mdot_str/4.0/M_PI/v_w_str;
const double k=0.0;
const double A=1.*M_PRO;
const double a=1.0; /* see Granot & Piran 2011 */
const double b=0.45; /* see Granot & Piran 2011 */

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

const int N_tbin = 2048;
const int N_ne = 2048;   /* need to be the same as N_nph */
const int N_nph = 2048; /* need to be the same as N_ne */

const int N_tobsbin = 128;
const double tobs_min=3.0e3;
const double tobs_max=3.0e8;

const double GAMMA_ELEC_MIN = 1.0;/* minimum Lorentz factor of electron */
const double GAMMA_ELEC_MAX = 2.0e9;/* maximum Lorentz factor of electron */
const double ENE_PH_MIN_DIV_MeC2 = 1.0e-14;/* mimimum energy of photon in unit of electron mass */
const double ENE_PH_MAX_DIV_MeC2 = 1.0e-4;/* maximum energy of photon in unit of electron mass */

/////////////////////////////////////////


#endif
