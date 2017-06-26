
/*******************************************************************************
SimFast21
Name: get_velocityfield                                                       
Calculates the velocity field from the linear density one
Outputs comoving velocity in Mpc (growth rate missing)
To get linear comoving peculiar velocity at a given z, just multiply by growth rate. You can use dgrowthdt(z) in auxiliary.c. The units will then be comoving Mpc/seconds.
*******************************************************************************/

/* --------------Includes ----------------------------------------- */
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <complex.h>   /* header for complex numbers in c */
#include <fftw3.h>     /* headers for FFTW library */

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){


  FILE *fid;
  size_t elem;
  long int i,j,p;
  long int indi, indj;
  double kk;
  fftwf_complex  *map_vel_c;
  float *map;
  fftwf_complex *map_in_c; 
  fftwf_plan pc2r;
  fftwf_plan pr2c;
  char fname[256];
  DIR* dir;

 
  if(argc != 2) {
    printf("Usage: get_velocityfield work_dir\n");
    printf("work_dir - directory containing simfast21.ini \n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);
  

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif

 if(!(map=(float *) fftwf_malloc(global_N_halo*global_N_halo*global_N_halo*sizeof(float)))) {  
    printf("Problem...\n");
    exit(1);
  }
 if(!(map_in_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
   printf(" Out of memory...\n");
   exit(1);
 }
 if(!(map_vel_c = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
   printf("Problem allocating memory for x velocity field in k-space...\n");
   exit(1);
 } 
 /* Fourier tranforms for the velocity boxes */
 if(!(pc2r=fftwf_plan_dft_c2r_3d(global_N_halo, global_N_halo, global_N_halo, map_vel_c, map, FFTW_ESTIMATE))) { 
   printf("Problem...\n");
   exit(1);
 }
 /* FFT to map */
 if(!(pr2c=fftwf_plan_dft_r2c_3d(global_N_halo, global_N_halo, global_N_halo, map , map_in_c, FFTW_ESTIMATE))) { 
   printf("Problem...\n");
   exit(1);
 }
 
 sprintf(fname, "%s/delta/delta_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,global_L/global_hubble);
 /* Read density file delta */
 fid=fopen(fname,"rb");	/* second argument contains name of input file */
 if (fid==NULL) {
   printf("\n Density file path is not correct or the file does not exit...\n"); 
   exit (1);
 }
 elem=fread(map,sizeof(float),global_N_halo*global_N_halo*global_N_halo,fid);
 fclose(fid);
 
 
 /***********************************************************************************/
 // FFT of density field (delta)
  
 fftwf_execute(pr2c);


 /********************************************************************/
 /********************************************************************/
 /********************************************************************/
 /* Computing velocity fields */

 /* Create directory Velocity */
  sprintf(fname,"%s/Velocity",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Velocity directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }
   
 /********************************************************************/
 printf("\nComputing v_x field...\n");fflush(0);

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N_halo,global_dk,map_vel_c,map_in_c,global_dx_halo) private(i,indi,j,indj,p,kk) 
#endif
 for(i=0;i<global_N_halo;i++) {         
   if(i>global_N_halo/2) {    /* Large frequencies are equivalent to smaller negative ones */
     indi=-(global_N_halo-i);
   }else indi=i;
   for(j=0;j<global_N_halo;j++) {           
     if(j>global_N_halo/2) {  
       indj=-(global_N_halo-j);
     }else indj=j;  
     for(p=0;p<=global_N_halo/2;p++) {
       kk=global_dk*sqrt(indi*indi+indj*indj+p*p);	
       if(kk>0){ 
	 // Normalize by including dx and dk
	 map_in_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=I*(global_dk)*(1/(kk*kk))*map_in_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]*global_dx_halo*global_dx_halo*global_dx_halo/global_L3;  
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=indi*map_in_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]/global_hubble; /* to turn velocity from Mpc/h to Mpc */  
       }else{
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=0;
       }
     }
   }
 }
 
 box_symmetriesf(map_vel_c,global_N_halo);
 
 /* Executes FFT */
 fftwf_execute(pc2r);
 
 printf("\nWriting v_x field to file...\n");fflush(0);

 /* velocity file in Mpc not Mpc/h */
 sprintf(fname, "%s/Velocity/vel_x_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,global_L/global_hubble); 
 if((fid=fopen(fname,"wb"))==NULL){  
   printf("\nThe file cannot be open\n");
   return 0;
 } 
 elem=fwrite(map,sizeof(float),global_N_halo*global_N_halo*global_N_halo,fid); 
 fclose(fid);


 /********************************************************************/
  printf("\nComputing v_y field...\n");fflush(0);

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N_halo,global_dk,map_vel_c,map_in_c,global_dx_halo) private(i,indi,j,indj,p,kk) 
#endif
 for(i=0;i<global_N_halo;i++) {         
   if(i>global_N_halo/2) {    /* Large frequencies are equivalent to smaller negative ones */
     indi=-(global_N_halo-i);
   }else indi=i;
   for(j=0;j<global_N_halo;j++) {           
     if(j>global_N_halo/2) {  
       indj=-(global_N_halo-j);
     }else indj=j;  
     for(p=0;p<=global_N_halo/2;p++) {
       kk=global_dk*sqrt(indi*indi+indj*indj+p*p);	
       if(kk>0){ 
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=indj*map_in_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]/global_hubble;  
       }else{
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=0;
       }
     }
   }
 }
 
 box_symmetriesf(map_vel_c,global_N_halo);
 
 /* Executes FFT */
 fftwf_execute(pc2r);
   
 printf("\nWriting v_y field to file...\n");fflush(0); 
 
 sprintf(fname, "%s/Velocity/vel_y_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,global_L/global_hubble); 
 if((fid=fopen(fname,"wb"))==NULL){  
   printf("\nThe file cannot be open\n");
   return 0;
 } 
 elem=fwrite(map,sizeof(float),global_N_halo*global_N_halo*global_N_halo,fid); 
 fclose(fid);


 /********************************************************************/
 printf("\nComputing v_z field...\n");fflush(0);

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N_halo,global_dk,map_vel_c,map_in_c,global_dx_halo) private(i,indi,j,indj,p,kk) 
#endif
 for(i=0;i<global_N_halo;i++) {         
   if(i>global_N_halo/2) {    /* Large frequencies are equivalent to smaller negative ones */
     indi=-(global_N_halo-i);
   }else indi=i;
   for(j=0;j<global_N_halo;j++) {           
     if(j>global_N_halo/2) {  
       indj=-(global_N_halo-j);
     }else indj=j;  
     for(p=0;p<=global_N_halo/2;p++) {
       kk=global_dk*sqrt(indi*indi+indj*indj+p*p);	
       if(kk>0){ 
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=p*map_in_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]/global_hubble;  
       }else{
	 map_vel_c[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=0;
       }
     }
   }
 }
 
 box_symmetriesf(map_vel_c,global_N_halo);
 
 /* Executes FFT */
 fftwf_execute(pc2r);
   
 printf("\nWriting v_z field to file...\n");fflush(0);
 
 sprintf(fname, "%s/Velocity/vel_z_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,global_L/global_hubble); 
 if((fid=fopen(fname,"wb"))==NULL){  
   printf("\nThe file cannot be open\n");
   return 0;
 } 
 elem=fwrite(map,sizeof(float),global_N_halo*global_N_halo*global_N_halo,fid); 
 fclose(fid);
 
 fftwf_free(map);
 fftwf_free(map_in_c);
 fftwf_free(map_vel_c);
 fftwf_destroy_plan(pc2r);
 fftwf_destroy_plan(pr2c);

 
 exit(0);    


}

