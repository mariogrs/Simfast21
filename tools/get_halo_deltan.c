
/*********************************************************************************************************
SimFast21
Auxiliar code - 2014
Description: Put the halo number fluctuations in a box for a given mass range. 
This is useful for instance to calculate the halo power spectrum and check it has the correct bias.
*********************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>  

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){

  
  FILE *fid;
  long int x,y,z;
  long int nhalos;
  float *halo_map;
  Halo_t *halo_v;
  size_t elem;
  char fname[300];
  double mass;
  long int index, i;
  double Mmin, Mmax, cat_Mmin, cat_Mmax, bin_Mmin, bin_Mmax; 
  double ntot;

  if(argc!=6) {
    printf("\nCalculates halo number density fluctuations in a box for a given mass range. Ouputs a box of NxNxN floats.\n");
    printf("usage: get_halo_deltan   Simfast21_dir   halo_catalog_file   Mmin   Mmax  output_box_file\n");
    printf("Halo catalog in Simfast21 format. Mmin and Mmax in Msun units\n");
    printf("Note: make sure that Mmin and Mmax is within the mass range of the halo finding algorithm!\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif

  
  Mmin=atof(argv[3]);
  Mmax=atof(argv[4]);
  sprintf(fname, "%s/Halos/%s",argv[1],argv[2]);
  
  if(!(halo_map=(float *) malloc(global_N3_halo*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  
    /* read halo catalog */
    if((fid=fopen(fname,"rb"))==NULL){  
      printf("Halo file: %s does not exist... Check path or run get_halos for this configuration\n",fname);
      exit(1);
    } 
    elem=fread(&nhalos,sizeof(long int),1,fid);
    printf("Reading %ld halos...\n",nhalos);fflush(0);
    if(!(halo_v=(Halo_t *) malloc(nhalos*sizeof(Halo_t)))) { 
      printf("Memory problem - halos...\n");
      exit(1);
    }
    elem=fread(halo_v,sizeof(Halo_t),nhalos,fid);  
    fclose(fid);
      
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(halo_map,global_N3_halo) private(i)
#endif
    for(i=0;i<global_N3_halo;i++)  halo_map[i]=0.0;

    printf("Searching through halo catalog...\n");fflush(0);
    ntot=0.0;
    cat_Mmin=(double)halo_v[0].Mass;
    cat_Mmax=(double)halo_v[0].Mass;
    bin_Mmin=cat_Mmin;
    bin_Mmax=cat_Mmax;
    for(i=0;i<nhalos;i++){
      mass=(double)halo_v[i].Mass;
      if(mass > cat_Mmax) cat_Mmax=mass;
      if(mass < cat_Mmin) cat_Mmin=mass;
      //      printf("halo mass: %E\n",mass);
      if(mass>=Mmin && mass <= Mmax) {	
	if(mass > bin_Mmax) bin_Mmax=mass;
	if(mass < bin_Mmin) bin_Mmin=mass;
	x= halo_v[i].x;     
	y= halo_v[i].y;  
	z= halo_v[i].z;    
	index=(long int)(x*global_N_halo*global_N_halo+y*global_N_halo+z);
	halo_map[index]+=1.0;
	ntot+=1.0;
      }
    }
    if((long int)ntot==0) {
      printf("No halos found in mass range. Exiting...\n");
      exit(1);
    }
    bin_Mmax=bin_Mmax*pow(10,global_halo_dlm); /* adjust mass interval due to mass resolution */
    printf("Catalogue Mmin: %E Msun, Catalogue Mmax: %E Msun\n",cat_Mmin, cat_Mmax);
    printf("Catalogue Mmin in given mass range: %E Msun, Catalogue Mmax in given mass range: %E Msun\n",bin_Mmin, bin_Mmax);
    printf("Total number of halos in mass range: %ld, average number of halos per cell: %E\n",(long int)ntot, ntot/global_N3_halo);
    printf("Number density: %E (h/Mpc)^3, dn/dm for given mass range: %E (h/Mpc)^3/Msun\n",ntot/global_L3,ntot/global_L3/(bin_Mmax-bin_Mmin)); fflush(0);

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(halo_map,global_N3_halo) private(i)
#endif
    for(i=0;i<global_N3_halo;i++)  halo_map[i]=halo_map[i]/(ntot/global_N3_halo)-1.0;
    
    printf("Writing file...\n"); fflush(0);
    sprintf(fname, "%s",argv[5]); 
    if((fid=fopen(fname,"wb"))==NULL){  
      printf("\nError opening output box\n");
      return 0;
    }
    elem=fwrite(halo_map,sizeof(float),global_N3_halo,fid);
    fclose(fid);
   
 exit(0);    
 
}
