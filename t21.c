
/*************************************************************
SimFast21
Calculates the 21cm brightness temperature
**NOTE**: output units in Kelvin
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


const double Tcmb0=2.725;       /* K */
const double maxcut=1.0;
const double mincut=-0.5;



int main(int argc, char * argv[]) {


  fftw_plan pc2r;
  fftw_plan pr2c;
  fftw_complex *dvdk;
  DIR* dir;
  float *temp;
  char fname[256];
  FILE *fid,*fidtxt, *fidtxt2;
  double *t21, *dvdr; 
  long int i,j,p,indi,indj;
  int nmaxcut,nmincut;
  double aver,Tcmb,kk,growth,dgdt,Hubz,TS,xtot;
  double zmin,zmax,dz,z;


  /* Check for correct number of parameters*/
  if(argc != 2) {
    printf("Usage : t21 base_dir\n");
    exit(1);
  }
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) printf("Lya and xray use set to false - assuming TS>>TCMB in t21 calculation.\n");
  zmin=global_Zminsim;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftw_init_threads();
  fftw_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
 

  /* Create directory deltaTb */
  sprintf(fname,"%s/deltaTb",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating t21 directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }

  sprintf(fname,"%s/Output_text_files/t21_av_N%ld_L%.1f.dat",argv[1],global_N_smooth,global_L/global_hubble); 
  if((fidtxt=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  
  sprintf(fname,"%s/Output_text_files/TS_av_N%ld_L%.1f.dat",argv[1],global_N_smooth,global_L/global_hubble); 
  if((fidtxt2=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  
  if(!(temp=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(t21=(double *) malloc(global_N3_smooth*sizeof(double)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(dvdk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*global_N_smooth*global_N_smooth*(global_N_smooth/2+1)))) {  
    printf("Problem...\n");
    exit(1);
  } 
  if(!(dvdr=(double *) fftw_malloc(global_N3_smooth*sizeof(double)))) {  
    printf("Problem...\n");
    exit(1);
  }
  if(!(pr2c=fftw_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, dvdr , dvdk, FFTW_MEASURE))) { 
    printf("Problem...\n");
    exit(1);
  }
  if(!(pc2r=fftw_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, dvdk, dvdr, FFTW_MEASURE))) { 
    printf("Problem...\n");
    exit(1);
  }

  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(z=zmax;z>(zmin-dz/10);z-=dz){    
    printf("T21 - z = %f\n",z);fflush(0);
    sprintf(fname,"%s/deltaTb/deltaTb_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
    if((fid = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping this redshift...\n",fname);
      fclose(fid);
    }else {      
      if(z>(global_Zminsfr-global_Dzsim/10) && global_use_Lya_xrays==1) {
	Tcmb=Tcmb0*(1.+z);
	sprintf(fname,"%s/x_c/xc_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
	if((fid = fopen(fname,"rb"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	fread(temp,sizeof(float),global_N3_smooth,fid);
	fclose(fid);
	for(i=0;i<global_N3_smooth;i++) t21[i]=(double)temp[i];
	sprintf(fname,"%s/Lya/xalpha_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
	if((fid = fopen(fname,"rb"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	fread(temp,sizeof(float),global_N3_smooth,fid);
	fclose(fid);
	aver=0.;
	for(i=0;i<global_N3_smooth;i++) {
	  xtot=(double)temp[i]+t21[i];
	  t21[i]=xtot/(1.+xtot);
	}
	
	sprintf(fname,"%s/xrays/TempX_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
	if((fid = fopen(fname,"rb"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}

	fread(temp,sizeof(float),global_N3_smooth,fid);
	fclose(fid);
	
	TS=0.;
	for(i=0;i<global_N3_smooth;i++) {
	
	  t21[i]=t21[i]*(1.-Tcmb/(double)temp[i]); // Temperature correction for high redshifts
	  TS+=Tcmb/(1.-t21[i]);
	}
      }else {
	for(i=0;i<global_N3_smooth;i++) t21[i]=1.0;
      }
      sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble); 
      fid=fopen(fname,"rb");
      if (fid==NULL) {printf("Error reading deltanl file... Check path or if the file exists..."); exit (1);}
      fread(temp,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
      for(i=0;i<global_N3_smooth;i++) dvdr[i]=(double)temp[i];
      /* Executes FFT */
      fftw_execute(pr2c);
      /********************************************************************/
      //     printf("Computing v field...\n"); 
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(dvdk,global_dk) private(i,indi,j,indj,p,kk) 
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
	      dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=0.;
	    }else {
	      kk=global_dk*sqrt(indi*indi+indj*indj+p*p);	
	      if(global_vi==1) {
		dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(global_dk*indi)*(global_dk*indi)/(kk*kk)*dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	      }else if(global_vi==2) {
		dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(global_dk*indj)*(global_dk*indj)/(kk*kk)*dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	      }else {
		dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=(global_dk*p)*(global_dk*p)/(kk*kk)*dvdk[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p];
	      }
	    }
	  }
	}
      }
      /* Executes FFT */
      fftw_execute(pc2r);
      nmaxcut=0;
      nmincut=0;
      growth=getGrowth(z);
      dgdt=dgrowthdt(z);
      Hubz=Hz(z);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(dvdr,global_L,global_dx_smooth,global_N3_smooth) private(i) 
#endif
      for(i=0;i<global_N3_smooth;i++) {
	dvdr[i]=-dgdt*dvdr[i]/global_L/global_L/global_L/Hubz*global_dx_smooth*global_dx_smooth*global_dx_smooth/growth; /* FT normalization + divide by H(z). Divide by growth to get deltanl back to z=0 */
      }
      for(i=0;i<global_N3_smooth;i++) {
	if(dvdr[i] > maxcut) {dvdr[i]=maxcut; nmaxcut++;}
	else if(dvdr[i] < mincut) {dvdr[i]=mincut; nmincut++;}
      }
      //      printf("nmaxcut: %d (%f %%), nmincut: %d (%f %%)\n\n",nmaxcut,100.*nmaxcut/global_N3_smooth,nmincut,100.*nmincut/global_N3_smooth);fflush(0);
     
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(t21,temp,dvdr,global_hubble,global_omega_b,global_omega_m,global_N3_smooth,z) private(i) 
#endif
      for(i=0;i<global_N3_smooth;i++) {
  	/* output in K */
	t21[i]=23.0/1000.*(1.+(double)temp[i])*t21[i]/(1.+1.*dvdr[i])*(0.7/global_hubble)*(global_omega_b*global_hubble*global_hubble/0.02)*sqrt((0.15/global_omega_m/global_hubble/global_hubble)*(1.+z)/10.);  /*Units in Kelvin */
      }
      
      sprintf(fname,"%s/Ionization/xHII_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((fid = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fread(temp,sizeof(float),global_N3_smooth,fid);  
      fclose(fid);
      for(i=0;i<global_N3_smooth;i++) {
	t21[i]=t21[i]*(1.-(double)temp[i]);  // neutral fraction...
      }
      
      sprintf(fname,"%s/deltaTb/deltaTb_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((fid = fopen(fname,"wb"))==NULL) {
	printf("Cannot open file:%s...\n",fname);
	exit(1);
      }
      for(i=0;i<global_N3_smooth;i++) temp[i]=(float)t21[i];
      fwrite(temp,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
      aver=0.;
      for(i=0; i<global_N3_smooth; i++) aver+=t21[i];
      fprintf(fidtxt,"%f %.8E\n",z,aver/global_N3_smooth);
      if(z>(global_Zminsfr-global_Dzsim/10) && global_use_Lya_xrays==1) fprintf(fidtxt2,"%f %.8E\n",z,TS/global_N3_smooth);
    }
  } /* ends z cycle */
    
  fclose(fidtxt);
  exit(0);

}


