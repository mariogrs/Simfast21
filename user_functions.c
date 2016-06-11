
/************************************************************************
SimFast21
User defined functions that are used in the simulation, such as fitting functions...
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

#define CM_PER_MPC 3.0857e24


/***************   Rion *************/
#define c1      7.25511524e+39
#define c2      9.367608e+07
#define c3      4.09039945e-01
#define c4      2.27044695e+00
#define V_norm  1.e+45

double Rion(float hmass, double redshift){
  double tmp_ion1 =  c1*hmass*pow((1.0 + redshift ),c4);
  double tmp_ion2 =  pow((hmass/c2),c3);
  double tmp_ion3 =  exp(- pow((c2/hmass),3.0));
  double tmp_ion4 =  tmp_ion1*tmp_ion2*tmp_ion3;
  double tmp_ion5 =  tmp_ion4*global_fesc/V_norm;
  return tmp_ion5;
}



/***************   Rrec *************/
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



/********** ratio between recombination rate coefficient (at T=10^4K) and interpolation function of Haardt & Madau (2012) uniform ionising background  ****************/

double G_H(double redshift){
  double gh_tmp1 =   1.86956756e-17*pow(redshift,5.)  - 9.05797228e-16*pow(redshift,4.) +  1.56916163e-14*pow(redshift,3.);
  double gh_tmp2 =   -1.07046180e-13*pow(redshift,2.) + 1.15923648e-13*redshift + 1.04351586e-12;
  double b_h =  4.19232273531e-13/(gh_tmp1 + gh_tmp2);
  return b_h;
}
/**************************** Popping et al. (2009) formula to compute the residual neutral fraction from ionising background ********/
double XHI(double ratio){
  double XHI_tmp1 = 2.*ratio + 1. - sqrt((2.*ratio + 1.)*(2.*ratio + 1.) - 4.*ratio*ratio  );
  double XHI_tmp3 = XHI_tmp1/(2.*ratio);
  if(ratio == 0.0) return 0.0; 
  else  return XHI_tmp3;
}

