
/*
SimFast21
Input_variables: Reads simfast21.ini
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "Input_variables.h"

                        			  
void get_Simfast21_params(char *basedir){
  
  FILE * paramfile;
  
  char line[256],fname[300];
  int length,n;
  char first[100], third[100];  
  char *cpoint;
 
  length=256;
  global_camb_file[0]=0;

  sprintf(fname, "%s/simfast21.ini",basedir); 
  if((paramfile=fopen(fname,"r"))==NULL){  
    printf("\nThe parameter file simfast21.ini cannot be open. Exit...\n");
    exit(1);
  }
  global_pk_flag=0;
  while(fgets(line,length,paramfile)!=NULL){
    if(line[0]!='#' && line[0]!='\n'){
      if((cpoint=strchr(line,'='))==NULL) {printf("Wrong format! Exiting...\n"); printf("%s\n",line); exit(1);}
      cpoint[0]=' ';
      n=sscanf(line,"%s %s",first, third);
      if(n!=2 || strlen(first)==0 || strlen(third)==0) {printf("Wrong format! Exiting...\n"); printf("%s\n",line); exit(1);}

      /* Simulation */
      else if(strcmp(first,"nthreads")==0)global_nthreads=atoi(third);
      else if(strcmp(first,"dz")==0)global_Dzsim=atof(third);
      else if(strcmp(first,"zmax")==0)global_Zmaxsim=atof(third);
      else if(strcmp(first,"zmin")==0) global_Zminsim=atof(third);
      else if(strcmp(first,"sim_length")==0)global_L=atof(third);  /* Mpc */
      else if(strcmp(first,"seed")==0){
	global_seed=atoi(third);
	if(global_seed < 1) {
	  printf("seed < 1 - no good: using time()\n");
	  global_seed=(long int)time(NULL);
	  printf("Using seed: %lu\n",global_seed);
	} 
      }     
      else if(strcmp(first,"Vel_comp")==0)global_vi=atoi(third);
      else if(strcmp(first,"N_halo")==0)global_N_halo=atoi(third);
      else if(strcmp(first,"N_smoothed")==0)global_N_smooth=atoi(third);

      /* Cosmology */
      if(strcmp(first,"use_camb_matterpower")==0){
	if(strcmp(third,"T")==0) global_pk_flag=1;    
	if(strcmp(third,"F")==0) global_pk_flag=0;    
      }
      else if(strcmp(first,"camb_file")==0)strcpy(global_camb_file,third);
      else if(strcmp(first,"omega_matter")==0)global_omega_m=atof(third);
      else if(strcmp(first,"omega_baryon")==0)global_omega_b=atof(third);
      else if(strcmp(first,"omega_lambda")==0)global_lambda=atof(third);
      else if(strcmp(first,"hubble")==0)global_hubble=atof(third);
      else if(strcmp(first,"spectral_index")==0)global_n_index=atof(third);
      else if(strcmp(first,"sigma8")==0)global_sig8_new=atof(third);    

      /* Halos */
      else if(strcmp(first,"critical_overdensity")==0)global_delta_c=atof(third);
      else if(strcmp(first,"STa")==0)global_STa=atof(third);
      else if(strcmp(first,"STb")==0)global_STb=atof(third);
      else if(strcmp(first,"STc")==0)global_STc=atof(third);
      else if(strcmp(first,"Use_subgrid")==0){
	if(strcmp(third,"T")==0) global_use_sgrid=1; else global_use_sgrid=0;
      }    
      else if(strcmp(first,"halo_Rmax")==0)global_halo_Rmax=atof(third);
      else if(strcmp(first,"halo_Rmin_dx")==0)global_halo_Rmin_dx=atof(third);
      else if(strcmp(first,"halo_Mmin")==0)global_halo_Mmin=atof(third);
      else if(strcmp(first,"halo_Nbins")==0)global_Nhbins=atoi(third);

      /* Ionization */
      else if(strcmp(first,"fesc")==0)global_fesc=atof(third);
      else if(strcmp(first,"Ion_cutoff")==0)global_xHlim=atof(third);
      else if(strcmp(first,"bubble_Rmax")==0)global_bubble_Rmax=atof(third);
      else if(strcmp(first,"bubble_Nbins")==0)global_bubble_Nbins=atoi(third);

      /* xray + Lya */
      else if(strcmp(first,"use_Lya_xrays")==0){
      	if(strcmp(third,"T")==0) global_use_Lya_xrays=1;    
	if(strcmp(third,"F")==0) global_use_Lya_xrays=0;
      }
      else if(strcmp(first,"Zminsfr")==0) global_Zminsfr=atof(third);
      else if(strcmp(first,"fstar")==0)global_fstar=atof(third);
      else if(strcmp(first,"Enu0")==0) global_Enu0=atof(third);
      else if(strcmp(first,"alpha_s")==0) global_alphas=atof(third);
      else if(strcmp(first,"L0")==0) global_L0=atof(third);
      else if(strcmp(first,"flux_Rmax")==0) global_flux_Rmax=atof(third); /* Mpc */
      else if(strcmp(first,"A_Lya")==0) global_A_Lya=atof(third);
      else if(strcmp(first,"alpha_Lya")==0) global_alpha_Lya=atof(third);

      /* Auxiliary */
      else if(strcmp(first,"Original_nldensity_box")==0){
	if(strcmp(third,"T")==0) global_save_original_deltanl=1;    
	if(strcmp(third,"F")==0) global_save_original_deltanl=0;
      }      
      else if(strcmp(first,"NL_halo_catalog")==0){
	if(strcmp(third,"T")==0) global_save_nl_halo_cat=1;    
	if(strcmp(third,"F")==0) global_save_nl_halo_cat=0;
      }

    }   
  }
  if(global_pk_flag==1) {
    printf("Using CAMB file for cosmology.\n");
    if(global_camb_file[0]==0) {
      printf("No CAMB file. Exiting...\n");
      exit(1);
    }
    sprintf(fname, "%s/%s",basedir,global_camb_file);
    printf("Using cosmology from CAMB file: %s\n",fname);
    set_cosmology_fromCAMB(fname);
  } else global_pk_flag=0;
  global_rho_m=global_omega_m*RHO_0/global_hubble;  /* units of M_sun/(Mpc/h)**3 */
  global_rho_b=global_omega_b*RHO_0/global_hubble;
  global_L=global_L*global_hubble;  /* Mpc/h */
  global_halo_Rmax=global_halo_Rmax*global_hubble;  /* Mpc/h */
  if(global_halo_Rmax>global_L/2.0) global_halo_Rmax=global_L/2.0;
  global_bubble_Rmax=global_bubble_Rmax*global_hubble;  /* Mpc/h */
  if(global_bubble_Rmax > global_L/2.) global_bubble_Rmax=global_L/2.; /* make sure maximum bubble radius is half the box size, e.g. one bubble over the all box... */
  global_flux_Rmax=global_flux_Rmax*global_hubble;  /* Mpc/h */
  global_L3=global_L*global_L*global_L;
  global_dk=2*PI/global_L;  /* h/Mpc */
  global_dx_halo=global_L/global_N_halo;
  global_dx_smooth=global_L/global_N_smooth;
  global_N3_halo=global_N_halo*global_N_halo*global_N_halo;
  global_N3_smooth=global_N_smooth*global_N_smooth*global_N_smooth;
  global_smooth_factor=global_N_halo/global_N_smooth;
  global_Zminsim=global_Zmaxsim-global_Dzsim*(double)ceil((global_Zmaxsim-global_Zminsim)/global_Dzsim); /* Make sure we include files for global_Zminsim */
  global_Zminsfr=global_Zmaxsim-global_Dzsim*(double)ceil((global_Zmaxsim-global_Zminsfr)/global_Dzsim);

  //  print_parms();

  fclose(paramfile);
}



void set_cosmology_fromCAMB(char * paramfilename) {

  FILE * paramfile;
  
  char line[256];
  int length,n;
  char first[99], third[99];  
  double ombh2,omch2,omk,w;
  char use_phys[99], root[99];
  char *cpoint;
  length=256;
 
  ombh2=0;
  omch2=0;
  omk=0;

  if((paramfile=fopen(paramfilename,"r"))==NULL){  
    printf("\nThe CAMB parameter file cannot be open. Exit...\n");
    exit(1);
  }
  while(fgets(line,length,paramfile)!=NULL){
    if(line[0]!='#' && line[0]!='\n'){
      if((cpoint=strchr(line,'='))==NULL) {printf("Wrong format! Exiting...\n"); printf("%s\n",line); exit(1);}
      cpoint[0]=' ';
      n=sscanf(line,"%s %s",first, third);
      if(n==0 || strlen(first)==0) {printf("Wrong format! Exiting...\n"); printf("%s\n",line); exit(1);}

      else if(strcmp(first,"use_physical")==0)strcpy(use_phys,third);
      else if(strcmp(first,"ombh2")==0)ombh2=atof(third);
      else if(strcmp(first,"omch2")==0)omch2=atof(third);
      else if(strcmp(first,"omk")==0)omk=atof(third);
      
      else if(strcmp(first,"omega_baryon")==0 && strcmp(use_phys,"F")==0)global_omega_b=atof(third);
      else if(strcmp(first,"omega_cdm")==0 && strcmp(use_phys,"F")==0)global_omega_m=atof(third);  /* CDM for now */
      else if(strcmp(first,"omega_lambda")==0 && strcmp(use_phys,"F")==0)global_lambda=atof(third);
      
      else if(strcmp(first,"hubble")==0)global_hubble=atof(third)/100.;
      else if(strcmp(first,"w")==0)w=atof(third);
      else if(strcmp(first,"scalar_spectral_index(1)")==0)global_n_index=atof(third);
      
      else if(strcmp(first,"output_root")==0)strcpy(root,third);
    
      else if(strcmp(first,"transfer_matterpower(1)")==0){
	strcat(root,"_");
	strcat(root,third);
	strcpy(global_pk_filename,root);
      }
    }
  }
  if(strcmp(use_phys,"T")==0){
    global_omega_b=ombh2/global_hubble/global_hubble;
    global_omega_m=omch2/global_hubble/global_hubble + global_omega_b;
    global_lambda=1.0-omk-global_omega_m;
  }  else global_omega_m=global_omega_m+global_omega_b;

}


void print_parms(void) {

  /*--------------------Simfast21 variables and parameters---------------------------*/
  printf("Input parameters:\n");
  printf("global_nthreads: %d\n",global_nthreads);
  printf("global_seed: %ld\n",global_seed);
  printf("global_N_halo: %ld\n",global_N_halo); //Linear number of cells of the box for determination of collapsed halos 
  printf("global_N3_halo: %ld\n",global_N3_halo); //Total number of cells of the box for determination of collapsed halos 
  printf("global_N_smooth: %ld\n",global_N_smooth); // Linear number of cells of the smoothed boxes  
  printf("global_N3_smooth: %ld\n",global_N3_smooth); // Total number of cells of the smoothed boxes  
  printf("global_smooth_factor: %f\n",global_smooth_factor); //Just N_halo/N_smooth
  printf("global_L: %f Mpc/h\n",global_L); //Physical size of the simulation box
  printf("global_L3: %f (Mpc/h)^3\n",global_L3);//Physical volume of the simulation box
  printf("global_dx_halo: %f Mpc/h\n",global_dx_halo);
  printf("global_dx_smooth: %f Mpc/h\n",global_dx_smooth);
  printf("global_dk: %f h/Mpc\n",global_dk);
  printf("global_vi: %d\n",global_vi);
  /* simulation range */
  printf("global_Dzsim: %f\n",global_Dzsim);
  printf("global_Zmaxsim: %f\n",global_Zmaxsim);
  printf("global_Zminsim: %f\n",global_Zminsim);
  printf("global_Zminsfr: %f\n",global_Zminsfr);

  /*-----------------------Cosmological parameters--------------------------- */
  if(global_pk_flag==0) printf("global_sig8_new: %f\n",global_sig8_new);
  printf("global_n_index: %f\n",global_n_index);
  printf("global_hubble: %f\n",global_hubble);
  printf("global_omega_m: %f\n",global_omega_m);
  printf("global_omega_b: %f\n", global_omega_b);  
  printf("global_lambda: %f\n",global_lambda);
  printf("global_rho_m: %E\n",global_rho_m);
  printf("global_rho_b: %E\n",global_rho_b);

  /*------------------------Halo collapse parameters-----------------------------*/
  printf("global_halo_Nbins: %d\n",global_Nhbins);
  printf("global_use_sgrid: %d\n",global_use_sgrid);
  printf("global_delta_c: %f\n",global_delta_c);
  printf("global_STa: %f\n",global_STa);
  printf("global_STb: %f\n",global_STb);
  printf("global_STc: %f\n",global_STc);
  printf("global_halo_Rmax: %f Mpc/h\n",global_halo_Rmax);
  printf("global_halo_Rmin_dx: %f (in cell units)\n",global_halo_Rmin_dx);
  printf("global_halo_Mmin: %E Msun\n",global_halo_Mmin);


  /*------------------------Ionization parameters-----------------------------*/
  printf("global_fesc: %f\n",global_fesc); //escape fraction
  printf("global_xHlim: %f\n",global_xHlim); //neutral fraction cutoff
  printf("global_bubble_Rmax: %f Mpc/h\n",global_bubble_Rmax);
  printf("global_bubble_Nbins: %d\n",global_bubble_Nbins);


  /*----------Variables for reading matter power spectrum from file-------- */
  printf("global_pk_flag: %d\n",global_pk_flag); // Matter power spectrum: 0 - Eisenstein & Hu fitting formulae; 1 - Read form file (Output of CMBFast, CAMB)
  printf("global_pk_filename: %s\n",global_pk_filename);
  printf("global_camb_file: %s\n",global_camb_file);


  /*-------------------Flags for output files and algorithm----------------------------------*/ 
  printf("global_save_nl_halo_cat: %d\n",global_save_nl_halo_cat);
  printf("global_save_original_deltanl: %d\n",global_save_original_deltanl);
  printf("global_use_Lya_xrays: %d\n",global_use_Lya_xrays);


  /*--------- Parameters for X-ray heating and Lya coupling----------*/
  printf("global_fstar: %f\n",global_fstar);
  printf("global_Enu0: %E\n",global_Enu0);
  printf("global_alphas: %f\n",global_alphas);
  printf("global_L0: %E\n",global_L0); 		
  printf("global_flux_Rmax: %f\n",global_flux_Rmax);
  printf("global_A_Lya: %f\n",global_A_Lya);
  printf("global_alpha_Lya: %E\n",global_alpha_Lya); 		
  printf("\n");

  fflush(0);

}

