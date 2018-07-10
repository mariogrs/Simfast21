
/*********************************************************************************************************
SimFast21
Auxiliar code - 2014
Description: Calculates halo dn/dm for a given halo catalogue.
Uses user provided logarithmic binning.
Also calculates theoretical mass function.
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
  double mass, dm,z, m1, m2, m3,dmi, dndma;
  long int ind, i, j;
  double Mmin, Mmax; 
  double dlm, smass, sdndma;
  double dndm[1000], massv[1000];
  int N;
  int nint=20;

  if(argc!=5) {
    printf("\nCalculates the halo dn/dm for a given catalogue.\n");
    printf("usage: get_dndm_nbins  Simfast21_dir   halo_catalog_file  z  nbins\n");
    printf("Halo catalog in Simfast21 format. Uses logarithmic binning.\n\n");
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
      
  printf("#Searching through halo catalog...\n");fflush(0);
  Mmin=(double)halo_v[0].Mass;
  Mmax=(double)halo_v[0].Mass;
  for(i=0;i<nhalos;i++){
    mass=(double)halo_v[i].Mass;
    if(mass > Mmax) Mmax=mass;
    if(mass < Mmin) Mmin=mass;
  }
  Mmax=Mmax*1.001;
  Mmin=Mmin*0.999;
  N=atoi(argv[4]);
  dlm=log10(Mmax/Mmin)/N;
  printf("# Mmax: %E    Mmin: %E     dlm: %E\n", Mmax, Mmin, dlm);
  printf("# N: %d\n", N);
  ntot=0;
  for(i=0;i<N;i++) {dndm[i]=0.0; massv[i]=0.0;}
  for(i=0;i<nhalos;i++){
    mass=(double)halo_v[i].Mass;
    ind=(int)((log10((mass)/Mmin)/dlm));
    dndm[ind]+=1.0;
    massv[ind]+=mass;
    ntot++;
  }
  if(ntot!=nhalos) {
    printf("Problem - halo numbers don't match!. Exiting...\n");
    exit(1);
  }
  printf("# Total number of halos in catalogue: %ld, average number of halos per cell: %E\n",ntot, 1.0*ntot/global_N3_halo);
  printf("# Number density: %E (h/Mpc)^3, dn/dm for total mass range: %E (h/Mpc)^3/Msun\n",1.0*ntot/global_L3,1.0*ntot/global_L3/(Mmax-Mmin));
  printf("\n# Mass [Msun]    dndm [(h/Mpc)^3/Msun]\n");
  printf("\n#  bin_min_Mass      bin_max_Mass      Mass_mid      Mass_av      Mass_weigth_sim      Mass_weight_calc      dndm_sim      dndm_calc\n");
  for(i=0;i<N;i++) {
    if(dndm[i]>0) massv[i]=massv[i]/dndm[i];
    m1=Mmin*pow(10,i*dlm);
    m2=Mmin*pow(10,(i+1)*dlm);
    m3=Mmin*pow(10,(i+1.0/2)*dlm);
    dmi=log(m2/m1)/nint;
    sdndma=0.0;
    smass=0;
    for(j=0;j<nint;j++) {
      mass=m1*exp(i*dmi+dmi/2.0);
      dndma=mass_function_ST(z,mass)*mass;  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume (times mass because of log integration) */
      sdndma+=dndma;  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume (times mass because of log integration) */
      smass+=mass*dndma;
    }
    smass=smass/sdndma;
    sdndma=sdndma*dmi;
    printf("%E     %E     %E     %E     %E     %E     %E     %E\n",m1, m2,m3, (m1+m2)/2.0, massv[i], smass, dndm[i]/global_L3/(m2-m1),sdndma/(m2-m1));
  }
  printf("\n");
  
  exit(0);}
  




