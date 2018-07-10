
/******************************************************************************************************************
Simple code to calculate analytical dndm...
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
  
  double z,dlm,m,dndm,sdndm, Mmin, Mmax;
  int i,N;

  if(argc!=6) {
    printf("\nCalculates the analytical dn/dm (halo mass function) for a given redshift and mass interval..\n");
    printf("usage: adndm   Simfast21_dir   z   Mmin   Mmax  N_bins\n");
    printf("Mmin and Mmax in Msun units. Uses logarithmic binning\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

  z=atof(argv[2]);
  Mmin=atof(argv[3]);
  Mmax=atof(argv[4]);
  N=atoi(argv[5]);

  dlm=log10(Mmax/Mmin)/N;
  sdndm=0.;
  printf("\n# Mass [Msun]     dndm [(h/Mpc)^3/Msun]\n");
  for(i=0;i<N;i++) {
    m=Mmin*pow(10,i*dlm+dlm/2.0);
    dndm=mass_function_ST(z,m);  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume */
    printf("%E          %E\n",m,dndm);
    sdndm+=dndm*m;
  }
  sdndm*=dlm;
  printf("halo number density over total mass range: %E (h/Mpc)^3\n",sdndm);
  printf("dndm over total mass range: %E (h/Mpc)^3/Msun\n",sdndm/(Mmax-Mmin));

}

