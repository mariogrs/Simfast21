
/******************************************************************************************************************
Simple code to calculate analytical z, r, nu
M Santos (2014)
*******************************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Input_variables.h"
#include "auxiliary.h"

const double nu0=1420.41;  /* MHz */

double get_r(double z);
double drdz(double z);


int main(int argc, char **argv){
  
  double z,zmin,zmax,dr,r,nu;

  if(argc!=5) {
    printf("\nOutputs: z, r (Mpc/h), nu (MHz)\n");
    printf("usage: rz.x   work_dir   zmin  zmax  dr\n");
    printf("work_dir is the folder with simfast21.ini. dr in Mpc/h\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

  zmin=atof(argv[2]);
  zmax=atof(argv[3]);
  dr=atof(argv[4]);
      
  z=zmin;
  for(r=get_r(zmin); r<=get_r(zmax); r+=dr) {
    nu=nu0/(1.0+z);
    printf("%f   %f   %f\n",z,r,nu);
    z=z+dr/drdz(z);
  }

}


/****************************************************************/

/* get r(z) in Mpc/h */
double get_r(double z) {
    
    double dz=0.001,r;
    int n,i;

    n=(int)(z/dz)+1;
    dz=z/n;

    r=0.;
    for(i=0; i<n;i++) {
	r+=drdz(i*dz+dz/2.);
    }

    return r*dz;
}


/* dr/dz in comoving Mpc/h */
double drdz(double z) {

    return 2997.9/sqrt(global_omega_m*(1.+z)*(1.+z)*(1.+z)+global_lambda);  /* value in Mpc/h */


}





