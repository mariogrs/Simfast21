
/*********************************************************************************************************
SimFast21
Auxiliar code - 2014
Description: Calculates halo dn/dm for a given halo catalogue. Uses same bins and mass range as simulation.
Also calculates theoretical mass function
*********************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>  

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){

  
  FILE *fid;
  long int nhalos, ntot;
  Halo_t *halo_v;
  size_t elem;
  char fname[300];
  double mass, dm,z, m1, m2, dmi, dndma,m3;
  long int ind, i, j;
  double Mmin, Mmax; 
  double dlm, R_lim, R, halo_mass;
  double dndm[10000];
  int N;
  int nint=20;

  if(argc!=4) {
    printf("\nCalculates the halo dn/dm for a given catalogue.\n");
    printf("usage: get_dndmb  work_dir   halo_catalog_file  z\n");
    printf("Halo catalog in Simfast21 format. Uses logarithmic binning and same mass bins as the simulation.\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);
  z=atof(argv[3]);
  sprintf(fname, "%s/Halos/%s",argv[1],argv[2]);
  
  /* read halo catalog */
  if((fid=fopen(fname,"rb"))==NULL){  
    printf("Halo file: %s does not exist... Check path or run get_halos for this configuration\n",fname);
    exit(1);
  } 
  elem=fread(&nhalos,sizeof(long int),1,fid);
  printf("#Reading %ld halos...\n",nhalos);fflush(0);
  if(!(halo_v=(Halo_t *) malloc(nhalos*sizeof(Halo_t)))) { 
    printf("Memory problem - halos...\n");
    exit(1);
  }
  elem=fread(halo_v,sizeof(Halo_t),nhalos,fid);  
  fclose(fid);

  R_lim=global_halo_Rmin_dx*global_dx_halo;
  R=R_lim*pow(global_Rhalo,(int)(log(global_halo_Rmax/R_lim)/log(global_Rhalo)))+R_lim/10000.;
  halo_mass=(4.0/3.0)*PI*global_rho_m*pow(R,3)*pow(global_Rhalo,3/2.0);  /* This is a new change ... */
  dlm=3*log10(global_Rhalo);
  N=(int)roundf((log10(halo_mass)-log10(global_halo_Mmin))/dlm)+1;
  Mmin=halo_mass/pow(global_Rhalo,3.0*N);
  if(Mmin >= global_halo_Mmin) {N=N+1; Mmin=Mmin/pow(global_Rhalo,3.0);}
  Mmin=Mmin/pow(global_Rhalo,3/2.0);

  printf("# Mmin: %E     dlm: %E\n", Mmin, dlm);
  printf("# N: %d\n", N);
  ntot=0;
  for(i=0;i<N;i++) dndm[i]=0.0;
  for(i=0;i<nhalos;i++){
    mass=(double)halo_v[i].Mass;
    ind=(int)roundf((log10((mass)/Mmin)/dlm));
    dndm[ind]+=1.0;
    ntot++;
  }
  if(ntot!=nhalos) {
    printf("Halo numbers don't match!. Exiting...\n");
    exit(1);
  }
  printf("# Total number of halos in catalogue: %ld, average number of halos per cell: %E\n",ntot, 1.0*ntot/global_N3_halo);
  printf("# Number density: %E (h/Mpc)^3, dn/dm for total mass range: %E (h/Mpc)^3/Msun\n",1.0*ntot/global_L3,1.0*ntot/global_L3/(Mmax-Mmin));
  printf("\n# Mass [Msun]    dndm [(h/Mpc)^3/Msun]\n");
  printf("\n#  Mass_1         Mass_2         Mass_3      dndm_sim       dndm_calc\n");
  for(i=0;i<N;i++) {
    m1=Mmin*pow(10,i*dlm);
    m2=Mmin*pow(10,(i+1)*dlm);
    m3=Mmin*pow(10,(i+1.0/2)*dlm);
    dmi=log(m2/m1)/nint;
    dndma=0.0;
    for(j=0;j<nint;j++) {
      mass=m1*exp(i*dmi+dmi/2.0);
      dndma+=mass_function_ST(z,mass)*mass;  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume (times mass because of log integration) */
    }
    dndma=dndma*dmi;
    printf("%E   %E   %E   %E  %E\n",m1, m2, m3, dndm[i]/global_L3/(m2-m1),dndma/(m2-m1));
  }
  printf("\n");
  
  exit(0);    

 }




