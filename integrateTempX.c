
/*
SimFast21
Calculates IGM temperature (output in Kelvin)
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
#define KK 1.3806580e-23   /* k in j/K */


double fh(double xi);
double fi(double xi);

double fdeltatime(double z1, double z2);

 
int main(int argc,char * argv[]) {

  float * xi, * xe, *TXold, *TX, *delta;
  double * Energy;
  FILE * filetx, *filetxold, *fid;
  FILE * filexi, *fileEnergycum, * fileEnergy,*fidtxt;
  char fname[256];
  double xitot;
  long int i,j;
  double z, hubble_factor;
  double deltatime;
  double zmin,zmax,dz,deltaTX,aver,zt,growth,growtht,dgdz;

  /* Check for correct number of parameters*/
  if(argc!=2) {
    printf("Integrates X ray temperature\n");
    printf("usage: integrateTempX base_dir\n");
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

  delta=(float*) malloc(sizeof(float)*global_N3_smooth);
  TX = (float*) malloc(sizeof(float)*Buffer);
  TXold = (float*) malloc(sizeof(float)*Buffer);
  xi = (float*) malloc(sizeof(float)*Buffer);
  xe = (float*) malloc(sizeof(float)*Buffer);
  Energy = (double*) malloc(sizeof(double)*Buffer);
  

  sprintf(fname,"%s/Output_text_files/Temp_av_N%ld_L%.1f.dat",argv[1],global_N_smooth,global_L/global_hubble); 
  if((fidtxt=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  
  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(z=zmax; z > (zmin-dz/10); z-=dz) {
    printf("z: %f\n",z);fflush(0);
    sprintf(fname,"%s/xrays/TempX_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
    if((filetx = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping redshift...\n",fname);
      fclose(filetx);
    }else {

      deltatime=fdeltatime(z,z+global_Dzsim);
      zt=z+global_Dzsim/2.;
      dgdz=dgrowthdz(zt);
      growth=getGrowth(z);
      growtht=getGrowth(zt);

      sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble); 
      fid=fopen(fname,"rb");
      if (fid==NULL) {printf("Error reading deltanl file... Check path or if the file exists..."); exit (1);}
      fread(delta,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
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
      if((fileEnergycum = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      sprintf(fname,"%s/xrays/TempX_z%.3f_N%ld_L%.1f.dat",argv[1],z,global_N_smooth,global_L/global_hubble);
      if((filetx = fopen(fname,"wb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      
      /* read previous temperature calculation at z+deltatime */    
      if (z<(zmax-dz/10)) {
	sprintf(fname,"%s/xrays/TempX_z%.3f_N%ld_L%.1f.dat",argv[1],z+global_Dzsim,global_N_smooth,global_L/global_hubble);    
	if((filetxold = fopen(fname,"rb"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	printf("Reading previous temp file\n");
      }else printf("no previous temp file\n");
      fflush(0);
      
      /* Assume adiabatic cooling between 1+z=100 and the highest redshift of the simulation (z+global_Dzsim) */
      /* Assume TK=180K at 1+z=100 (decoupled from CMB) */
      /* Note: if starting redshift is too low, adiabatic cooling assumption may be wrong */ 
      /* This factor includes the integration in 2/3dn/dt/n */
      hubble_factor=1+2.*Hz(z)*deltatime;
      aver=0.;
      for (i=0;i<Nbufs;i++) {
	fread(xi,sizeof(float),Buffer,filexi);  
	fread(xe,sizeof(float),Buffer,fileEnergycum); /* ionization due to xrays */
        fread(Energy,sizeof(double),Buffer,fileEnergy);  /* epsilonX over n */
	if (z<(zmax-dz/10)) fread(TXold,sizeof(float),Buffer,filetxold);
	else {
	  for (j=0;j<Buffer;j++){ 
            TXold[j]=180.*pow((1.+z+global_Dzsim)/100.,2);  /* Kelvin */
	    TXold[j]=TXold[j]*pow(1+delta[i*Buffer+j],2./3);
          }
	}
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(TX,deltatime,TXold,hubble_factor,xe,xi,Energy,global_Dzsim,zt,delta,z,growth,growtht,dgdz) private(j,deltaTX,xitot) 
#endif
	for (j=0;j<Buffer;j++) {
	  /* Add a minimum ionization state of the gas
	     due to X-ray to prevent total non-heating
	     from X-ray.
	  */
	  xitot= xe[j] + xi[j];
	  if (xitot > 1.) xitot=1.;
	  if (xitot < 0.) xitot=0.;
	  deltaTX=Energy[j]*fh(xitot)*2./3./KK;  /* DTX in Kelvin/s */
	  if(isnan(deltaTX)) {
	    printf("Error nan: %ld %E %E %E %E %E\n",j,Energy[j],fh(xitot),xitot,xe[j],xi[j]);
	    exit(1);
	  }
	  TX[j]=(float)((TXold[j]*(1-global_Dzsim/(1.+zt)-global_Dzsim*dgdz*delta[i*Buffer+j]/3.0/(growth+growtht*delta[i*Buffer+j]))+deltaTX*deltatime)/(1+global_Dzsim/(1.+zt)+global_Dzsim*dgdz*delta[i*Buffer+j]/3.0/(growth+growtht*delta[i*Buffer+j])));
        }
	for (j=0;j<Buffer;j++) aver+=TX[j];
	fwrite(TX,sizeof(float),Buffer,filetx);
      }
      fprintf(fidtxt,"%f %.8E\n",z,aver/global_N3_smooth);
      fclose(filetx);
      fclose(filexi);
      fclose(fileEnergycum);
      fclose(fileEnergy);
      if (z<(zmax-dz/10)) fclose(filetxold);
    }
  } /* ends z cycle */  
  fclose(fidtxt);
  free(TX);
  free(TXold);
  free(xi);
  free(xe);
  free(Energy);
    
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


