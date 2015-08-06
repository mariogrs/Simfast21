
/******************************************************************************************************************
Simple code to calculate the analytical halo bias...
M Santos (2014)
*******************************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Input_variables.h"
#include "auxiliary.h"



int main(int argc, char **argv){
  
  double z,dlm,m,sdndm, Mmin, Mmax;
  int i,N;
  double temp, bias, R, sigmaM;

  if(argc!=5) {
    printf("\nCalculates the analytical halo bias for a given redshift and mass interval.\n");
    printf("usage: abias   work_dir   z   Mmin   Mmax\n");
    printf("Mmin and Mmax in Msun units.\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

  z=atof(argv[2]);
  Mmin=atof(argv[3]);
  Mmax=atof(argv[4]);

  dlm=0.01;
  N=log10(Mmax/Mmin)/dlm;
  sdndm=0.;
  bias=0.0;
  for(i=0;i<N;i++) {
    m=Mmin*pow(10,i*dlm+dlm/2.0);
    R=pow(m*3.0/4.0/PI/global_rho_m,1./3);
    sigmaM=(double)sigma(R);
    temp=mass_function_ST(z,m)*m; /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume */
    sdndm+=temp;
    bias=bias+temp*Bias(z,sigmaM);
  }
  sdndm*=dlm;
  bias=bias*dlm/sdndm;
  printf("halo number density over mass range: %E (h/Mpc)^3\n",sdndm);
  printf("bias over mass range: %f\n",bias);

}

