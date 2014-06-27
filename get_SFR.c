
/*********************************************************************************************************
SimFast21
Description: Calculates SFR from the collapsed mass field
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



int main(int argc, char **argv){

  
  FILE *fid;
  DIR* dir;
  char fname[300];
  long int i,j,nz,elem;
  double aver[1000],fcollv[1000],sfrv[1000],zv[1000];
  float *halo_mass;
  double zmin,zmax,dz,redshift;
  
  if(argc==1 || argc > 5) {
    printf("Generates SFRD using nonlinear halo boxes\n");
    printf("usage: get_SFR base_dir [zmin] [zmax] [dz]\n");
    printf("base_dir contains simfast21.ini\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) {printf("Lya and xray use set to false - no need to calculate SFR\n");exit(0);}
  if(argc > 2) {
    zmin=atof(argv[2]);
    if (zmin < global_Zminsfr) zmin=global_Zminsfr;
    if(argc > 3) {
      zmax=atof(argv[3]);
      if(zmax>global_Zmaxsim) zmax=global_Zmaxsim;
      if(argc==5) dz=atof(argv[4]); else dz=global_Dzsim;
    }else {
      zmax=global_Zmaxsim;
      dz=global_Dzsim;
    }
    zmin=zmax-dz*ceil((zmax-zmin)/dz); /* make sure (zmax-zmin)/dz is an integer so that we get same redshifts starting from zmin or zmax...*/ 
  }else {
    zmin=global_Zminsfr;
    zmax=global_Zmaxsim;
    dz=global_Dzsim;
  }
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

  if(!(halo_mass=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }

  nz=0;
  printf("Calculating halo mass average...\n");
  for(redshift=zmin;redshift<(zmax+dz/10);redshift+=dz){
    zv[nz]=redshift;
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
    if((fid=fopen(fname,"rb"))==NULL){  
      printf("\nError opening %s\n",fname);
      exit(1);
    }
    elem=fread(halo_mass,sizeof(float),global_N3_smooth,fid);
    fclose(fid);
    aver[nz]=0.;
    for(i=0;i<global_N3_smooth;i++) aver[nz]+=halo_mass[i];
    fcollv[nz]=aver[nz]/(global_rho_m*global_L3);  /* Msun/(Mpc/h)^3 */
    aver[nz]/=global_N3_smooth;
    nz++;
  }  

  printf("Interpolation for SFRD...\n");
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nz);
  gsl_spline_init (spline, zv, fcollv, nz);


  sprintf(fname,"%s/Output_text_files",argv[1]);
  if((dir=opendir(fname))==NULL) {
    printf("Creating Output_text_files directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }
  sprintf(fname,"%s/Output_text_files/sfrd_av_N%ld_L%.0f.dat",argv[1],global_N_smooth,global_L); 
  if((fid=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  
  for(i=nz-1;i>=0;i--) {
    sfrv[i]=gsl_spline_eval_deriv(spline, zv[i], acc)*dzdt(zv[i])/3.171E-8*global_rho_b*global_fstar; /* SFRD in Msun/(Mpc/h)^3/year */
    if(sfrv[i]<0.) sfrv[i]=0.;
    fprintf(fid,"%f %.8E\n",zv[i],sfrv[i]);
  }  
  fclose(fid);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  printf("Writing SFRD files...\n");
  for(i=0;i<nz;i++) {
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.0f.dat",argv[1],zv[i],global_N_smooth,global_L); 
    if((fid=fopen(fname,"rb"))==NULL){  
      printf("\nError opening %s\n",fname);
      exit(1);
    }
    elem=fread(halo_mass,sizeof(float),global_N3_smooth,fid);
    fclose(fid);
    
    if(aver[i]>0.)
      for(j=0;j<global_N3_smooth;j++) halo_mass[j]*=sfrv[i]/aver[i];
    else for(j=0;j<global_N3_smooth;j++) halo_mass[j]=0.;
    sprintf(fname, "%s/SFR/sfrd_z%.3f_N%ld_L%.0f.dat",argv[1],zv[i],global_N_smooth,global_L); 
    if((fid=fopen(fname,"wb"))==NULL){
      printf("\nError opening output sfrd file... Chech if path is correct...\n"); 
    }
    elem=fwrite(halo_mass,sizeof(float),global_N3_smooth,fid);
    fclose(fid); 
  }

  exit(0);

}  
  
