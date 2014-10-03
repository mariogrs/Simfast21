
/*********************************************************************************************************
SimFast21
Auxiliar code - 2014
Description: Calculates halo dn/dm for a given halo catalogue
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
  long int nhalos;
  Halo_t *halo_v;
  size_t elem;
  char fname[300];
  double mass, dm;
  long int ind, i;
  double Mmin, Mmax; 
  double ntot,sim_Mmax,sim_Mmin, dlm;
  double dndm[1000];
  int N;

  if(argc!=3) {
    printf("\nCalculates the halo dn/dm for a given catalogue.\n");
    printf("usage: get_dndm   work_dir   halo_catalog_file\n");
    printf("Halo catalog in Simfast21 format. Uses logarithmic binning and same mass bins as the simulation.\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

  sprintf(fname, "%s/Halos/%s",argv[1],argv[2]);
  sim_Mmax=(4.0/3.0)*PI*global_rho_m*pow(global_halo_Rmax,3);
  sim_Mmin=(4.0/3.0)*PI*global_rho_m*pow(global_halo_Rmin_dx*global_dx_halo,3);
  //  printf("#Input: %s, output: %s, Mmin: %E Msun, Mmax: %E Msun, Sim Mmin: %E Msun, Sim Mmax: %E Msun\n",fname, argv[5], Mmin, Mmax, sim_Mmin, sim_Mmax);

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
    dlm=3*log10(global_Rhalo);
    Mmax=Mmax*pow(global_Rhalo,3);
    N=(log10(Mmax)-log10(Mmin))/dlm+1;
    ntot=0.0;
    for(i=0;i<N;i++) dndm[i]=0.0;
    for(i=0;i<nhalos;i++){
      mass=(double)halo_v[i].Mass;
      ind=(int)roundf(log10(mass/Mmin)/dlm);
      dndm[ind]+=1.0;
      ntot+=1.0;
    }
    if((long int)ntot==0) {
      printf("No halos found in mass range. Exiting...\n");
      exit(1);
    }
    printf("# Total number of halos in mass range: %ld, average number of halos per cell: %E\n",(long int)ntot, ntot/global_N3_halo);
    printf("# Number density: %E (h/Mpc)^3, dn/dm for given mass range: %E (h/Mpc)^3/Msun\n",ntot/global_L3,ntot/global_L3/(Mmax-Mmin));
    printf("\n# Mass [Msun]    dndm [(h/Mpc)^3/Msun]\n");
    for(i=0;i<N;i++) {
      mass=Mmin*pow(10,i*dlm+dlm/2.0);
      dm=Mmin*pow(10,i*dlm)*(pow(10,dlm)-1.0);
      printf("%E      %E\n",mass, dndm[i]/global_L3/dm);
    }
    printf("\n");
   
 exit(0);    
 
}
