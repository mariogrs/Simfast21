
/* 
SimFast21
Integrates eq. 24 in Santos et al 2008 paper
eg calculates xe (ionization only due to X rays
also uses ionization from bubbles to get fion 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _OMPTHREAD_
#include <omp.h>
#endif   
#include <gsl/gsl_integration.h>

#include "Input_variables.h"
#include "auxiliary.h"

#define Buffer global_N_smooth*global_N_smooth
#define Nbufs global_N_smooth

#define C  0.9971
#define A 0.2663
#define B  1.3163

#define C2 0.3908
#define A2 0.4092
#define B2 1.7592


#define EV2J 1.60217646e-19


double fh(double xi);
double fi(double xi);

double fdeltatime(double z1, double z2);

 
int main(int argc,char * argv[]) {

  float *xi, *Energycumio, *Energycumioold;
  double *Energy, deltatime,Energy2io;
  FILE *filexi, *fileEnergycum, *fileEnergycumold, *fileEnergy;
  long int i,j;
  char fname[1000];
  double zmin,zmax,dz,z,xitot;


  /* Check for correct number of parameters*/
  if(argc!=2) {
    printf("Integrates xe for X ray temperature calculation\n");
    printf("usage: integratexe base_dir\n");
    printf("base_dir contains simfast21.ini\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) {printf("Lya and xray use set to false - no need to calculate xrays\n");exit(0);}
  zmin=global_Zminsfr;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
  

  xi = (float*) malloc(sizeof(float)*Buffer);
  Energycumio = (float*) malloc(sizeof(float)*Buffer);
  Energy = (double*) malloc(sizeof(double)*Buffer);  /* note the double */
  Energycumioold = (float*) malloc(sizeof(float)*Buffer);
  Energy2io=13.6*EV2J;    
  
  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(z=zmax; z >(zmin-dz/10); z-=dz) {
    printf("z: %f\n",z);fflush(0);
    sprintf(fname,"%s/xrays/xe_heat_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
    if((fileEnergycum = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping this redshift...\n",fname);
      fclose(fileEnergycum);
    }else {
      deltatime=fdeltatime(z,z+global_Dzsim);
      //      printf("deltatime : %E seconds\n",deltatime);    
      sprintf(fname,"%s/Ionization/xHII_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((filexi = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      sprintf(fname,"%s/xrays/EpsilonXon_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((fileEnergy = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      sprintf(fname,"%s/xrays/xe_heat_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((fileEnergycum = fopen(fname,"wb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      if (z<(zmax-dz/10)) {
	sprintf(fname,"%s/xrays/xe_heat_z%.3f_N%ld_L%.1f.dat",argv[1],z+global_Dzsim,global_N_smooth,global_L/global_hubble);
	if((fileEnergycumold = fopen(fname,"rb"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	printf("Reading previous xe file\n");
      } else printf("no previous xe file\n");
      fflush(0);
      
      memset(xi,0,sizeof(float)*Buffer);
      memset(Energy,0,sizeof(double)*Buffer);
      memset(Energycumioold,0,sizeof(float)*Buffer);
      
      printf("Looping over box for z=%f\n",z);fflush(0);
      for (i=0;i<Nbufs;i++) {   
	/* xi from bubble box */
	fread(xi,sizeof(float),Buffer,filexi);
	fread(Energy,sizeof(double),Buffer,fileEnergy);  /* epsilonX over n - in double */
	if (z<(zmax-dz/10)) 
	  fread(Energycumioold,sizeof(float),Buffer,fileEnergycumold); /* read previous xe */
	/* Assume xe=0 for larger z */
	
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Energy,Energycumio,Energycumioold,xi,deltatime,Energy2io) private(j,xitot) 
#endif
	for (j=0;j<Buffer;j++) {
	  /* Add a minimum ionization state of the gas
	     due to X-ray to prevent total non-heating
	     from X-ray.
	  */
	  xitot=xi[j]+Energycumioold[j];  /* adding previous xe to make it a better approximation */
	  if(xitot>1.0) xitot=1.0;      
	  if(xitot<0.) xitot=0.0;      
	  Energycumio[j]=(float)(Energy[j]*deltatime/Energy2io*fi(xitot)+ Energycumioold[j]);
	}
	
	fwrite(Energycumio,sizeof(float),Buffer,fileEnergycum);
      }
      if (z<(zmax-dz/10)) fclose(fileEnergycumold);
      fclose(fileEnergy);
      fclose(filexi);
      fclose(fileEnergycum);
    }
  } /* ends z cycle */    
    free(xi);
    free(Energycumio);
    free(Energy);
    free(Energycumioold);
    
    exit(0);
  }




double dtdz(double z) {
  return (double)(1./dzdt(z));
}


double ddtdz(double z, void * params) {
  return dtdz(z);
}



double fdeltatime(double z1, double z2) {

  double temp;
  double dz=0.001,z;
  double deltat;
  int n,i;

  if (z1 > z2) {
    temp=z1;
    z1=z2;
    z2=temp;
  }
  n=(int)round((z2-z1)/dz);
  dz=(z2-z1)/n;
  deltat=0.;
  for (i=0;i<n;i++) {
    z=z1+dz/2.+dz*i;
    deltat+=fabs(dtdz(z));  // just in case there is a negative sign  -  results in seconds
  }
  deltat*=dz;

  return deltat;
}



double fh (double xi) {
    return C*(1-pow(1-pow(xi,A),B));
}

double fi (double xi) {
    return C2*pow(1-pow(xi,A2),B2);
}


