
/*
SimFast21
Calculates optical depth from the ionization files
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Input_variables.h"

double sigt=6.652e-25;  /* cm^2 */


/* dr/dz in comoving Mpc/h */
double drdz(double z) {

    return 2997.9/sqrt(global_omega_m*(1.+z)*(1.+z)*(1.+z)+global_lambda);  /* value in Mpc/h */


}


int main(int argc, char * argv[]) {

  int i,nz,ind;
  double z,dz=0.01,zv[1000],xH[1000];
  double nH0,nHe0,ne0,tau,tau2,ddz;
  FILE *fid;
  char fname[300];

  /* Check for correct number of parameters*/
  if (argc != 2) {
    printf("Usage : opdepth base_dir (for simulation)\n");
    exit(0);
  }

  get_Simfast21_params(argv[1]);

  sprintf(fname,"%s/Output_text_files/zsim.txt",argv[1]);
  if((fid = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  i=0;
  while(fscanf(fid,"%lf",&(zv[i]))==1) i++;
  nz=i;
  printf("Number of lines: %d\n",nz);
  fclose(fid);    
  sprintf(fname, "%s/Output_text_files/x_HI_N%ld_L%.1f.dat",argv[1],global_N_smooth,global_L/global_hubble);
  if((fid = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  i=0;
  while(fscanf(fid,"%lf",&(xH[i]))==1) i++;
  fclose(fid);
  if(i!=nz) {printf("Error!\n");exit(1);}
  nH0=(1.-YHe)*(global_rho_b*Msun/pow(Mpc2m*100./global_hubble,3))/mH; /* cm^-3 */
  nHe0=(YHe)*(global_rho_b*Msun/pow(Mpc2m*100./global_hubble,3))/mHe;  /* cm^-3 */
  ne0=nH0+nHe0;
  tau=0.;
  for(z=0.;z<zv[nz-1];z+=dz) {
    tau+=(1.+z)*(1.+z)*drdz(z);
  }
  tau=tau*dz*ne0*sigt*Mpc2m/global_hubble*100.;
  printf("Contribution after universe is fully ionised: %E\n",tau);
  ddz=(zv[0]-zv[nz-1])/(nz-1);
  tau2=0.;
  while(z<=zv[0]) {
    ind=(int)round((z-zv[nz-1])/ddz);
    if(ind>nz-1 || ind <0) printf("Error - ind: %d\n",ind);
    ind=nz-1-ind;
    tau2+=(1.+z)*(1.+z)*drdz(z)*(1.-xH[ind]);
    z+=dz;
  }
  tau2=tau2*dz*ne0*sigt*Mpc2m/global_hubble*100.;
  printf("Contribution from the EoR: %E\n",tau2);
  printf("Total: %E\n",tau+tau2);
  
  exit(0);

}



