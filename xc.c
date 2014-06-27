/*************************************************************
SimFast21
Calculates xc (collisional coupling)
*************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Input_variables.h"


double kappa_HI(double T);
double kappa_e(double T);

const double A10=2.85e-15;      /* s^-1 */
const double Tstar=0.068;       /* K */ 
const double Tcmb0=2.725;       /* K */



int main(int argc, char * argv[]) {


  DIR* dir;
  double nH_bar0, nHe_bar0,nH_bar, nHe_bar,xc_HI,xc_e,aver,Tcmb;
  float *xc, *temp, *fHII, *density_map;
  char fname[256];
  FILE * fid,*fidtxt;
  double * zbox; 
  long int i;
  int nzsfr,nzbox;
  double xHII;

  /* Check for correct number of parameters*/
  if (argc != 2) {
    printf("Usage : xc base_dir\n");
    exit(0);
  }
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) {printf("Lya and xray use set to false - no need to calculate xc\n");exit(0);}

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
 

  /* Create directory xc */
  sprintf(fname,"%s/x_c",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating xc directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }

  nzsfr=(int)((global_Zmaxsim-global_Zminsfr)/global_Dzsim)+1;
  zbox = (double*) malloc(sizeof(double)*nzsfr);

  sprintf(fname,"%s/Output_text_files/zsim.txt",argv[1]);
  if((fid = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  for (i=0;i<nzsfr;i++) 
    if(fscanf(fid,"%lf",&(zbox[nzsfr-1-i]))!=1) { 
      printf("Error reading zsim.txt!\n");
      exit(1);
    }
  fclose(fid);    
  sprintf(fname,"%s/Output_text_files/xc_av_N%ld_L%.0f.dat",argv[1],global_N_smooth,global_L); 
  if((fidtxt=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  
  if(!(density_map=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(fHII=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(temp=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  if(!(xc=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  nH_bar0=(1.-YHe)*(global_rho_b*Msun/pow(Mpc2m*100./global_hubble,3))/mH; /* cm^-3 */
  nHe_bar0=(YHe)*(global_rho_b*Msun/pow(Mpc2m*100./global_hubble,3))/mHe;  /* cm^-3 */
  printf("nH0: %.6E  nHe0: %.6E\n",nH_bar0,nHe_bar0);

  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(nzbox=nzsfr-1;nzbox >= 0;nzbox--) {

    /* our box is at nzbox */
    printf("ztocompute: %f\n",zbox[nzbox]);fflush(0);
    sprintf(fname,"%s/x_c/xc_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
    if((fid = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping this redshift...\n",fname);
      fclose(fid);
    }else {
      sprintf(fname,"%s/Ionization/xHII_z%.3f_eff%.2lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_eff,global_N_smooth,global_L);
      if((fid = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fread(fHII,sizeof(float),global_N3_smooth,fid);  
      fclose(fid);
      sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L); 
      fid=fopen(fname,"rb");
      if (fid==NULL) {printf("\nError reading deltanl file... Check path or if the file exists..."); exit (1);}
      fread(density_map,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
      sprintf(fname,"%s/xrays/TempX_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
      if((fid = fopen(fname,"rb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fread(temp,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
      
      nH_bar=nH_bar0*pow(1.+zbox[nzbox],3);
      nHe_bar=nHe_bar0*pow(1.+zbox[nzbox],3);         
      Tcmb=Tcmb0*(1.+zbox[nzbox]);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(xc,global_N3_smooth,fHII,density_map,temp,Tcmb,nH_bar,nHe_bar) private(i,xc_HI,xc_e) 
#endif
      for(i=0; i<global_N3_smooth; i++) {
	xHII=(double)fHII[i];
	if(xHII < 0.99999) 
	  xc_HI=4./3*(double)(1.-xHII)*nH_bar*(1.+(double)density_map[i])*kappa_HI(temp[i])/A10*Tstar/Tcmb;
	else xc_HI=0.;
	if(xHII>0.00001) 
	  xc_e=4./3*((double)xHII*nH_bar+nHe_bar*(double)xHII)*(1.+(double)density_map[i])*kappa_e(temp[i])/A10*Tstar/Tcmb;
	else xc_e=0.;
	//    xc_e=((1-(double)fHI[i])*nH_bar+nHe_bar*(2.-2.*fHeI[i]-fHeII[i]))*(double)rho_mat[i]*kappa_e(temp)/A10*Tstar/Tcmb;
	xc[i]=(float)(xc_HI+xc_e);
      }
      sprintf(fname,"%s/x_c/xc_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
      if((fid = fopen(fname,"wb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fwrite(xc,sizeof(float),global_N3_smooth,fid);
      fclose(fid);
      aver=0.;
      for(i=0; i<global_N3_smooth; i++) aver+=xc[i];
      fprintf(fidtxt,"%f %.6E\n",zbox[nzbox],aver/global_N3_smooth);
    }
  } /* ends z cycle */
    
  fclose(fidtxt);
  exit(0);

}





double kappa_HI(double T) {
  /* T in Kelvin */

  return (3.1e-11*pow(T,0.357)*exp(-32./T));  /* cm^3 s^-1 */

}

double kappa_e(double T) {
    /* T in Kelvin */

    double logT;

    //    if(T>1.0e4) T=1.0e4;
    logT=log10(T);
    if(logT < 0.) logT=0.;

    return pow(10.,(-9.607+0.5*logT*exp(-pow(logT,4.5)/1800.)));  /* cm^3 s^-1 */

}

