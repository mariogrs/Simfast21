
/*************************************************************
Routines to get velocities and gradients...
Outputs: 
  v/c (unitless - proper peculiar velocity/speed of light)
  (dv/ds)/H (unitless - proper gradient of the proper peculiar velocity along the given direction, divided by H(z))
*************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <complex.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fftw3.h>     /* headers for FFTW library */

#include "Input_variables.h"
#include "auxiliary.h"

int main(int argc, char * argv[]) {


  fftw_plan pc2r1;
  fftw_plan pc2r2;
  fftw_plan pr2c1;
  fftw_complex *dvdk, *velk;
  float *temp;
  char fname[256];
  FILE *fid;
  double *vel, *dvdr; 
  long int i,j,p,indi,indj;
  int vi;
  double kk,growth,dgdt,Hubz,z;
  double c=299792458.0; /* m/s */

  /* Check for correct number of parameters*/
  if(argc != 4) {
    printf("Usage : vel_grad base_dir z axis\nbase_dir: directory where Simfast21 simulation is.\nz: redshift where to calculate the velocities.\naxis: direction along which to calculate the gradient (1- X; 2-Y; 3-Z) - Z is fastest moving index in boxes.\n");
    printf("Outputs:\n  v/c (unitless - proper peculiar velocity/speed of light)\n  (dv/ds)/H (unitless - proper gradient of the proper peculiar velocity along the given direction, divided by H(z))\n");
    exit(1);
  }
  get_Simfast21_params(argv[1]);
  z=atof(argv[2]);
  vi=atoi(argv[3]);

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftw_init_threads();
  fftw_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
 

  if(!(temp=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(velk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*global_N_smooth*global_N_smooth*(global_N_smooth/2+1)))) {  
    printf("Problem...\n");
    exit(1);
  } 
  if(!(dvdk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*global_N_smooth*global_N_smooth*(global_N_smooth/2+1)))) {  
    printf("Problem...\n");
    exit(1);
  } 
  if(!(vel=(double *) fftw_malloc(global_N3_smooth*sizeof(double)))) {  
    printf("Problem...\n");
    exit(1);
  }
  if(!(dvdr=(double *) fftw_malloc(global_N3_smooth*sizeof(double)))) {  
    printf("Problem...\n");
    exit(1);
  }
  if(!(pr2c1=fftw_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, vel , velk, FFTW_MEASURE))) { 
    printf("Problem...\n");
    exit(1);
  }
  if(!(pc2r1=fftw_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, velk, vel, FFTW_MEASURE))) { 
    printf("Problem...\n");
    exit(1);
  }
  if(!(pc2r2=fftw_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, dvdk, dvdr, FFTW_MEASURE))) { 
    printf("Problem...\n");
    exit(1);
  }

  sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble); 
  fid=fopen(fname,"rb");
  if (fid==NULL) {printf("Error reading deltanl file... Check path or if the file exists..."); exit (1);}
  fread(temp,sizeof(float),global_N3_smooth,fid);   /* read density field */
  fclose(fid);
  for(i=0;i<global_N3_smooth;i++) vel[i]=(double)temp[i];
  /* Executes FFT */
  fftw_execute(pr2c1);   /* FT density field */
  /********************************************************************/
  printf("Computing v field...\n"); 
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(dvdk,global_dk,velk) private(i,indi,j,indj,p,kk) 
#endif
  for(i=0;i<global_N_smooth;i++) {         
    if(i>global_N_smooth/2) {    /* Large frequencies are equivalent to smaller negative ones */
      indi=-(global_N_smooth-i);
    }else indi=i;
    for(j=0;j<global_N_smooth;j++) {           
      if(j>global_N_smooth/2) {  
	indj=-(global_N_smooth-j);
      }else indj=j;  
      for(p=0;p<=global_N_smooth/2;p++) {
	if(i==0 && j==0 && p==0) {
	  velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=0.;
	  dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=0.;
	}else {
	  kk=global_dk*sqrt(indi*indi+indj*indj+p*p);	
	  if(vi==1) {
	    velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*indi)/(kk*kk)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	    dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*indi)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	  }else if(vi==2) {
	    velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*indj)/(kk*kk)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	    dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*indj)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	  }else {
	    velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*p)/(kk*kk)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	    dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(I*global_dk*p)*velk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	  }
	}
      }
    }
  }
  /* Executes FFT */
  fftw_execute(pc2r1);
  fftw_execute(pc2r2);
  growth=getGrowth(z);
  dgdt=dgrowthdt(z);
  Hubz=Hz(z);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(dvdr,global_L,global_dx_smooth,global_N3_smooth,vel) private(i) 
#endif
  for(i=0;i<global_N3_smooth;i++) {
    vel[i]=dgdt*vel[i]/global_L/global_L/global_L*global_dx_smooth*global_dx_smooth*global_dx_smooth/growth/(1+z)/c*3.08567758e22/global_hubble; /* FT normalization /(1+z) to get proper, last bit to get units of m/s (3.08567758e22 to get Mpc in m and 1/h because velocity is in Mpc/h initially). Divide by growth to get deltanl back to z=0 */
    dvdr[i]=dgdt*dvdr[i]/global_L/global_L/global_L/Hubz*global_dx_smooth*global_dx_smooth*global_dx_smooth/growth; /* FT normalization + divide by H(z). Divide by growth to get deltanl back to z=0 */
  }

  /* This outputs the velocity in proper units (hence the divide by (1+z)), divided by c: v/c in draft */
  sprintf(fname, "%s/Velocity/vel_z%.3f_N%ld_L%.1f_%d.dat",argv[1],z,global_N_smooth,global_L/global_hubble,vi);
  if((fid = fopen(fname,"wb"))==NULL) {
    printf("Cannot open file:%s...\n",fname);
    exit(1);
  }
  printf("Outputing %s\n",fname);
  for(i=0;i<global_N3_smooth;i++) temp[i]=(float)vel[i];
  fwrite(temp,sizeof(float),global_N3_smooth,fid);
  fclose(fid);
  
  /* This outputs dvc/dr (all comoving) or dv/ds (all proper - same thing) divided by H: dv/ds/H in the draft */
  sprintf(fname, "%s/Velocity/dvdr_z%.3f_N%ld_L%.1f_%d.dat",argv[1],z,global_N_smooth,global_L/global_hubble,vi);
  if((fid = fopen(fname,"wb"))==NULL) {
    printf("Cannot open file:%s...\n",fname);
    exit(1);
  }
  printf("Outputing %s\n",fname);
  for(i=0;i<global_N3_smooth;i++) temp[i]=(float)dvdr[i];
  fwrite(temp,sizeof(float),global_N3_smooth,fid);
  fclose(fid);
    
  exit(0);

}


