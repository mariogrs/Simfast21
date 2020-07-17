
/************************************************************************
SimFast21
User defined functions that are used in the simulation, such as fitting functions...
See arXiv:1510.04280 for more details
*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>   /* header for complex numbers in c */

#include "Input_variables.h"
#include "auxiliary.h"

double Rion(float hmass, double redshift);
double Rrec(float overdensity, double redshift);
double G_H(double redshift);
double XHI(double ratio);
double sfr(float hmass, double z);
double Qion(double z);
double ebG(double nu);


#define CM_PER_MPC 3.0857e24


/* units of s^-1 */
/***************   Rion *************/
/* ionisation rate */
//#define c1      7.25511524e+39
#define c1      1.575e39
#define c2      9.367608e+07
//#define c3      4.09039945e-01
#define c3      0.44
#define c4      2.27044695e+00
#define V_norm  1.e+45

double Rion(float hmass, double redshift){
  double tmp_ion1 =  c1*hmass*pow((1.0 + redshift ),c4);
  double tmp_ion2 =  pow((hmass/c2),c3);
  double tmp_ion3 =  exp(- pow((c2/hmass),3.0));
  double tmp_ion4 =  tmp_ion1*tmp_ion2*tmp_ion3;
  double tmp_ion5 =  tmp_ion4/V_norm;
  return tmp_ion5;
}



/***************   Rrec *************/
/* recombination rate */
#define r1 9.84815696
#define r2 1.76105133
#define r3 0.81722878
#define r4 5.06773134

double Rrec(float overdensity, double redshift){
  double rec1 =  1.e-24*r1*pow((1.0 + redshift ),r4);
  double rec2 =  pow((overdensity/r2),r3);
  double rec3 =  pow(rec2/(1.+rec2),4.);
  double DVV = pow(global_L*CM_PER_MPC/global_N_smooth/global_hubble/(1. + redshift),3.0);
  double rec55 = rec1*rec3*DVV/V_norm;
  return rec55;
}


/********* SFR ****************/
/* units: M/yr */

double sfr(float hmass, double z) {

  return (Rion(hmass,z)/Qion(z));

}


/********* Qion ****************/
/* units: yr*s^-1*M^-1 */

#define a 2.05707274e-02
#define b 5.29944413e+01

double Qion(double z) {

  return pow(10.0, a*z + b)/V_norm;

}


/********** ratio between recombination rate coefficient (at T=10^4K) and interpolation function of Haardt & Madau (2012) uniform ionising background  ****************/
double G_H(double redshift){
  double gh_tmp1 =   - 4.37152818e-05*pow(redshift,5.)  + 1.64394635e-03*pow(redshift,4.) - 2.03044756e-02*pow(redshift,3.);
  double gh_tmp2 =   7.18335221e-02*pow(redshift,2.) - 8.85791388e-02*redshift -1.20887500e+01;
  double b_h =  4.19232273531e-13/pow(10.0, gh_tmp1 + gh_tmp2);
  return b_h;
}


/**************************** Popping et al. (2009) formula to compute the residual neutral fraction from ionising background ********/
double XHI(double ratio){
  double XHI_tmp1 = 2.*ratio + 1. - sqrt((2.*ratio + 1.)*(2.*ratio + 1.) - 4.*ratio*ratio  );
  double XHI_tmp3 = XHI_tmp1/(2.*ratio);
  if(ratio == 0.0) return 0.0; 
  else  return XHI_tmp3;
}



/* for xalpha calculation
/* generic function for number photons/baryon/freq.(Hz)  - use: A*nu^-alpha */
/* Assumes A/nu^0.9 SED with 20,000 Lyman photons/baryon between Lya and Ly_limit */
/* We can change A and index 0.9 (A is degenerate with SFR efficiency */
/* Frequency in Hz */
double ebG(double nu) {

  double A_Lya = 1979.878526;
  double alpha_Lya = 0.9;

  return A_Lya*pow(nu,-alpha_Lya);

}
