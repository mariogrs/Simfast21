
/*********************************************************************************************************
SimFast21
Description: Calculates SFR from the non-linear halo field using a fitting function calibrated to simulations
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
     
#include "Input_variables.h"
#include "auxiliary.h"

double sfr(float hmass, double z);


int main(int argc, char **argv){

  
  FILE *fid, *fid_sa;
  DIR* dir;
  char fname[300];
  long int i,j,nz;
  double sfrd_aver;
  float *sfr_box1, *sfr_box2;
  double zmin,zmax,dz,redshift;
  long int nhalos;
  Halo_t *halo_v;
  size_t elem;

  
  if(argc !=2) {
    printf("Generates SFRD using nonlinear halo boxes\n");
    printf("usage: get_SFR base_dir\n");
    printf("base_dir contains simfast21.ini\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  if((global_use_Lya_xrays==0) && (global_use_SFR==0)) {printf("Lya and xray use set to false - no need to calculate SFR\n");exit(0);}
  if(global_use_SFR==1) { /* use SFR all the way... */
    zmin=global_Zminsim;
  } else zmin=global_Zminsfr;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;
  
  printf("\nCalculating SFRD between z=%f and z=%f with step %f\n",zmin,zmax,dz);
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
  
  /* Create directory SFR */
  sprintf(fname,"%s/SFR",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating SFR directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }  
  if(!(sfr_box1=(float *) malloc(global_N3_halo*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(sfr_box2=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }

  sprintf(fname,"%s/Output_text_files",argv[1]);
  if((dir=opendir(fname))==NULL) {
    printf("Creating Output_text_files directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }
  sprintf(fname,"%s/Output_text_files/sfrd_av_N%ld_L%.1f.dat",argv[1],global_N_smooth,global_L/global_hubble); 
  if((fid_sa=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  

  
  /* redshift cycle */
  for(redshift=zmin;redshift<(zmax+dz/10);redshift+=dz){
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,global_L/global_hubble); 
    fid=fopen(fname,"rb");
    if (fid==NULL) {printf("\nError reading %s file... Check path or if the file exists...",fname); exit (1);}
    elem=fread(&nhalos,sizeof(long int),1,fid);
    printf("Reading %ld halos...\n",nhalos);fflush(0);
    if(!(halo_v=(Halo_t *) malloc(nhalos*sizeof(Halo_t)))) { 
      printf("Problem - halo...\n");
      exit(1);
    }
    elem=fread(halo_v,sizeof(Halo_t),nhalos,fid);
    fclose(fid);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_halo, sfr_box1) private(i)
#endif
    for(i=0;i<(global_N3_halo);i++){
      sfr_box1[i] =0.0;
    }     
    // CIC smooth Rion//
    for(i=0;i<nhalos;i++){
      CIC_smoothing(halo_v[i].x, halo_v[i].y, halo_v[i].z, sfr(halo_v[i].Mass, redshift), sfr_box1, global_N_halo); /* SFR in Msun/yr */
    }
    free(halo_v);
    smooth_sum(sfr_box1, sfr_box2, global_N_halo, global_N_smooth); /* add SFR over smoothed cells */
    /* convert to SFRD in units of Msun/(Mpc/h)^3/year */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_smooth, sfr_box2) private(i)
#endif
    for(i=0;i<(global_N3_smooth);i++){
      sfr_box2[i] = sfr_box2[i]/global_dx_smooth/global_dx_smooth/global_dx_smooth;
    }     
    printf("Writing SFRD file...\n");
    sprintf(fname, "%s/SFR/sfrd_z%.3f_N%ld_L%.1f.dat",argv[1],redshift,global_N_smooth,global_L/global_hubble);
    fid=fopen(fname,"wb");
    if (fid==NULL) {printf("\nError opening %s file...\n",fname); exit (1);}
    elem=fwrite(sfr_box2,sizeof(float),global_N3_smooth,fid);
    fclose(fid); 

  /* average */
  sfrd_aver=0.0;
  for(i=0;i<(global_N3_smooth);i++){
    sfrd_aver+=sfr_box2[i];
  }
  sfrd_aver=sfrd_aver/global_N3_smooth;
  fprintf(fid_sa,"%f %.8E\n",redshift,sfrd_aver); /* prints average SFRD */

  }  /* ends redshift cycle */
  
  fclose(fid_sa); 
  exit(0);

}  
  
