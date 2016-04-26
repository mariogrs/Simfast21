
/*********************************************************************************************************
SimFast21
Auxiliar code - 2014
Description: Halo format conversion...
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
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){

  FILE *fid;
  DIR* dir;
  int nhalos;
  long int ln;
  int ix,iy,iz,i;
  Halo_sim *halo_in;
  Halo_t *halo_out;
  char fname[300];
  size_t elem;
  double z, dx;

  if(argc!=4) {
    printf("\nConverts halo catalogues from Mellema to Simfast21 format\n");
    printf("usage: convert_halos.x  work_dir  halo_catalog_input  z\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);
  print_parms();

  z=atof(argv[3]);
  sprintf(fname, "%s/Halos/halo_z%.3f_N%ld_L%.0f.dat.catalog",argv[1],z,global_N_halo,global_L);
  fid=fopen(fname,"rb");
    if(fid!=NULL) {
      printf("File:%s already exists - exiting...\n",fname); 
      exit(1);
    }

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif

  dx=global_dx_halo/global_hubble; /* convert to Mpc */

  /* read halo catalog */
  if((fid=fopen(argv[2],"rb"))==NULL){  
    printf("Halo file: %s does not exist...\n",argv[2]);
    exit(1);
  } 
  elem=fread(&nhalos,sizeof(int),1,fid);
  printf("Reading %d halos...\n",nhalos);fflush(0);
  if(!(halo_in=(Halo_sim *) malloc(nhalos*sizeof(Halo_sim)))) { 
    printf("Memory problem - halos...\n");
    exit(1);
  }
  elem=fread(halo_in,sizeof(Halo_sim),nhalos,fid);  
  if((int)elem!=nhalos) {printf("Problem reading halos. Exiting...\n"); exit(1);}
  fclose(fid);
  if(!(halo_out=(Halo_t *) malloc(nhalos*sizeof(Halo_t)))) { 
    printf("Memory problem - halos...\n");
    exit(1);
  }
   
  printf("Converting...\n"); fflush(0);
  //#ifdef _OMPTHREAD_
  //#pragma omp parallel for shared(halo_in,halo_out) private(i)
  //#endif
  for(i=0;i<nhalos;i++)  {
    halo_out[i].Mass=pow(10.0,halo_in[i].val[10]);  /* Mass in Msun */
    halo_out[i].Radius=halo_in[i].val[9]*(1.0+z)/1000.0*global_hubble;  /* Convert from radius in proper Kpc to comoving Mpc/h*/
    ix=(int)roundf(halo_in[i].val[0]/dx);
    if(ix < 0 || ix >= global_N_halo) {
      //#pragma omp critical
      {
	printf("Error in halo index: i, ix = %d, %d. Exiting...\n",i,ix);
	exit(1);
      }
    }
    iy=(int)roundf(halo_in[i].val[1]/dx);
    if(iy < 0 || iy >= global_N_halo) {
      //#pragma omp critical
      {
	printf("Error in halo index: i, iy = %d, %d. Exiting...\n",i, iy);
	exit(1);
      }
    }
    iz=(int)roundf(halo_in[i].val[2]/dx);
    if(iz < 0 || iz >= global_N_halo) {
      //#pragma omp critical
      {
	printf("Error in halo index: i, iz = %d, %d. Exiting...\n",i, iz);
	exit(1);
      }
    }
    /* switch indexes since we're assuming x should be the slowest moving index, not the fastest as in Garrelt... */
    halo_out[i].x=iz;
    halo_out[i].y=iy;
    halo_out[i].z=ix;
  }
  
  printf("Writing...\n"); fflush(0);
  sprintf(fname,"%s/Halos",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Halo directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }  
  sprintf(fname, "%s/Halos/halo_z%.3f_N%ld_L%.0f.dat.catalog",argv[1],z,global_N_halo,global_L);
  fid=fopen(fname,"wb");
  if(fid==NULL){printf("\nCatalog file output error. Exiting...\n"); exit (1);}
  ln=(long int)nhalos;
  elem=fwrite(&ln,sizeof(long int),1,fid);
  elem=fwrite(halo_out,sizeof(Halo_t),nhalos,fid);
  if((int)elem!=nhalos) {printf("Problem writing halos. Exiting...\n"); exit(1);}
  
  exit(0);    
  
}
