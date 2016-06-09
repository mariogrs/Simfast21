
/******************************************************************************************************************
SimFast21
Desciption: This routines determines collapsed halos positions and masses
from the linear density field. The barrier condition is given by Seth & Thoren modified version
of the critical overdensity.
Halos with mass smaller than cell size determined using a subgrid method
******************************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>   /* header for complex numbers in c */
#include <fftw3.h>     /* headers for FFTW library */
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Input_variables.h"
#include "auxiliary.h"



int main(int argc, char **argv){

  
  double dm, mass_aux, dndm_av, hbias, tempv;
  double rhalo;
  long int nhalos_z, nhalos_R;
  double Pnhalos;
  int nhalos_cat, nn;
  fftwf_plan pc2r;
  fftwf_plan pr2c;
  fftwf_complex *map_in;
  fftwf_complex *map_in_aux;
  FILE *fid_in;
  FILE *fid_out;
  size_t elem;
  double redshift;
  long int i,j,p,ii,jj,pp,indi,indj;
  long int a,b,c;
  int flag;
  long int ncells1D;
  long int ind;
  long int index, indice;
  long int *halo_ind;
  char  *flag_halo;
  double deltaCMZ;
  double halo_mass;
  double R_lim;
  Halo_t *halo_cat;
  float *map_in_f, *mass;
  double sigma_aux;
  double kk;
  double R2,z2,x2,y2;
  char halo_filename[300];
  char dens_filename[300];
  char fname[256];
  DIR* dir;
  double zmax,zmin,dz;
  double growth,R;
  gsl_rng * rpoisson = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set (rpoisson, 0);
  
  if(argc != 2) {
    printf("Generates Halo catalogs for a range of redshifts\n");
    printf("Usage: get_halos work_dir\n");
    printf("work_dir - directory containing simfast21.ini \n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  //  print_parms();
  
  zmin=global_Zminsim;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;
    
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);fflush(0);
#endif

  if(global_halo_Rmin_dx < 2) {
    printf("Warning: \"halo_Rmin_dx\" is quite small - this might have resolution effects on the standard halo finding excursion set formalism\n");
    R_lim=1.0*global_dx_halo;
    //    R_lim=0.620350491*global_dx_halo;
  }else {
    R_lim=global_halo_Rmin_dx*global_dx_halo;
  }
  halo_mass=(4.0/3.0)*PI*global_rho_m*pow(R_lim,3);
  printf("Minimum mass for standard halo finding method (excursion set formalism): %E\n",halo_mass);
  if(global_use_sgrid==1){
    if(halo_mass <= global_halo_Mmin+10.) {
      printf("No need to do subgridding: resolution is enough for %E mass halos\n",global_halo_Mmin);
      global_use_sgrid=0;
      R_lim=pow(3./4/PI/global_rho_m*global_halo_Mmin,1./3);
    }
  }else {
    if(halo_mass > global_halo_Mmin) {
      printf("Warning - minimum mass for simulation is larger than halo_Mmin: you need to use subgrid.\n");
    } else R_lim=pow(3./4/PI/global_rho_m*global_halo_Mmin,1./3);
  }
  if(global_use_sgrid==1)
    R=pow(3./4/PI/global_rho_m*global_halo_Mmin,1./3);
  else
    R=R_lim;
  rhalo=exp(log(global_halo_Rmax/R)/global_Nhbins);
  //  printf("%f  %f  %f  %f  %f  %d\n",global_halo_Rmax, global_L,R, R_lim, rhalo, global_Nhbins); fflush(0);
  
  if(!(halo_cat=(Halo_t *) fftwf_malloc(global_N_halo*global_N_halo*sizeof(Halo_t)))) {    /* get memory for float map */
    printf("Memory problem: halo_cat\n");
   exit(1);
  } 
  if(!(map_in_f=(float *) fftwf_malloc(global_N3_halo*sizeof(float)))) {    /* get memory for float map */
    printf("Problem1...\n");
   exit(1);
  } 
  if(!(mass=(float *) fftwf_malloc(global_N3_halo*sizeof(float)))) {    /* get memory for float map */
    printf("Problem1...\n");
   exit(1);
  } 
  if(!(map_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
    printf("Problem2...\n");
    exit(1);
  }
  if(!(pr2c=fftwf_plan_dft_r2c_3d(global_N_halo, global_N_halo, global_N_halo, map_in_f , map_in, FFTW_ESTIMATE))) { 
    printf("Problem3...\n");
    exit(1);
  }    

  /* Reading dark matter density box */
  printf("Reading dark matter density box\n");fflush(0);
  sprintf(dens_filename, "%s/delta/delta_z0_N%ld_L%.1f.dat", argv[1],global_N_halo, global_L/global_hubble);  
  fid_in=fopen(dens_filename,"rb");	/* second argument contains name of input file */
  if (fid_in==NULL) {
    printf("Error reading density file %s - check if the file exists...\n",dens_filename); 
    exit (1);
  }  
  elem=fread(map_in_f,sizeof(float),global_N3_halo,fid_in);    
  fclose(fid_in);
  fftwf_execute(pr2c);
  
  if(!(map_in_aux = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
    printf("Problem5...\n");
    exit(1);
  }
  if(!(flag_halo=(char *) malloc(global_N3_halo*sizeof(char)))) {  
    printf("Problem6...\n");
    exit(1);
  }
  if(!(halo_ind=(long int *) malloc(global_N3_halo/200*sizeof(long int)))) { 
    printf("Problem7...\n");
    exit(1);
  }
  
  //  printf("smooth_f: %ld\n",smooth_factor);  
  //  printf("Preparing FFT from filtered k-space to real space...\n");fflush(0);
  if(!(pc2r=fftwf_plan_dft_c2r_3d(global_N_halo, global_N_halo, global_N_halo,(fftwf_complex *) map_in_aux,(float *)map_in_aux , FFTW_MEASURE ))) { 
    printf("Problem...\n");
    exit(1);
  }
  //  printf("Done...\n");fflush(0);
  sprintf(fname,"%s/Halos",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Halo directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }  

  /****************************************************/
  /********  Redshift cycle ********/
  printf("Redshift cycle...\n");fflush(0);
  for(redshift=zmax;redshift>(zmin-dz/10);redshift-=dz){
    
    printf("\n\n\n +++++++++++++++++++++++++++ z = %f +++++++++++++++++++++++++++++++\n",redshift);fflush(0);
    sprintf(halo_filename, "%s/Halos/halo_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,(global_L/global_hubble));
    fid_out=fopen(halo_filename,"rb");
    if(fid_out!=NULL) {
      printf("File:%s already exists - skipping this redshift\n",halo_filename); fflush(0);
      fclose(fid_out);
    }else {

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(mass, flag_halo,global_N_halo) private(i)
#endif
      for(i=0;i<global_N3_halo;i++)  {flag_halo[i]=0; mass[i]=0.0;}
  
    growth=getGrowth(redshift);        
    sprintf(halo_filename, "%s/Halos/halo_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,(global_L/global_hubble));
    fid_out=fopen(halo_filename,"wb");
    if(fid_out==NULL){printf("\n Catalog file error\n"); exit (1);}
    nhalos_z=0;  /* counts total halos for one box */
    elem=fwrite(&nhalos_z,sizeof(long int),1,fid_out); /* reserve space for number of halos at beginning of file */    
    R=global_halo_Rmax;
    halo_mass=(4.0/3.0)*PI*global_rho_m*(pow(R,3)+pow(R/rhalo,3))/2.0;  /* This is a new change - it seems this mass fits theory better... */
 
    /**************** start halo cycle for one box ****************/
    while(R>=R_lim && halo_mass>=global_halo_Mmin){   

      printf("\n\n----------------------------------------- R=%f  M=%E ----------------------------------------\n",R,halo_mass);fflush(0);
      R2=R*R;   
      sigma_aux=sigma(R);
      deltaCMZ=growth*deltaFilter(sigma_aux,growth);
      if ((sigma_aux*growth*6) < deltaCMZ){    /* just in case probability is very low */
	printf("sigma(R=%f) << delta_c=%f, skipping \n",R,deltaCMZ);      
      } else {
	
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_in_aux,map_in,global_N_halo,global_dx_halo,R) private(i,j,p,indi,indj,kk)
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
	      /* 3d vector k (frequency) is just (indi, indj, p)*global_dk (global_dk is the moduli of the unit vector) */
	      kk=global_dk*sqrt(indi*indi+indj*indj+p*p); 
	      //map_in_aux[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]=map_in[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]*W_filter(kk*R)*dx*dx*dx;
	      *((fftwf_complex *)map_in_aux + (i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p))= *((fftwf_complex *)map_in + (i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p))*W_filter(kk*R)*global_dx_halo*global_dx_halo*global_dx_halo/global_L3*growth;  /* take care of all norms now */	      
	    }
	  }
	}	
	printf("Executing FFT to real space...\n");fflush(0);
      /* Executes FFT */
	fftwf_execute(pc2r);  /* Inline transform - careful with output format... */
	
	nhalos_R=0; /* halo count for one size R */	
	ncells1D=(long int)(R/global_dx_halo);
	if(ncells1D==0)ncells1D=1;      
	printf("Checking for halos and overlap...\n");fflush(0);
	for(i=0;i<global_N_halo;i++){ 
	  nhalos_cat=0;
	  for(j=0;j<global_N_halo;j++){    
	    for(p=0;p<global_N_halo;p++){
	      index=i*2*global_N_halo*(global_N_halo/2+1)+j*2*(global_N_halo/2+1)+p;  	  
	      if(*((float *)map_in_aux + index)>deltaCMZ){  /* found a halo... */
		flag=0;
		ind=0;
		for(ii=-ncells1D;ii<=ncells1D && flag==0;ii++){
		  x2=ii*ii*global_dx_halo*global_dx_halo;
		  a=i+ii;
		  a=check_borders(a,global_N_halo);
		  for(jj=-ncells1D;jj<=ncells1D && flag==0;jj++){
		    y2=jj*jj*global_dx_halo*global_dx_halo; 
		    if(x2+y2<=R2){		
		      b=j+jj;
		      b=check_borders(b,global_N_halo);
		      for(pp=-ncells1D;pp<=ncells1D && flag==0;pp++){
			z2=pp*global_dx_halo*pp*global_dx_halo;	  
			if(x2+y2+z2<=R2){
			  c=p+pp;
			  c=check_borders(c,global_N_halo);
			  if(flag_halo[a*global_N_halo*global_N_halo+b*global_N_halo+c]!=1){
			    halo_ind[ind]=a*global_N_halo*global_N_halo+b*global_N_halo+c;
			    ind++;
			  }else flag=1;
			}
		      }
		    }				
		  }
		}
		if (flag==0){
		  nhalos_R++;
		  nhalos_z++;
		  for(indice=0;indice<ind;indice++)flag_halo[halo_ind[indice]]=1;
		  halo_cat[nhalos_cat].x=i;
		  halo_cat[nhalos_cat].y=j;
		  halo_cat[nhalos_cat].z=p;
		  halo_cat[nhalos_cat].Mass=halo_mass; /* mass in Msun */
		  //		  halo_cat[nhalos_cat].Mass=halo_mass*(1+*((float *)map_in_aux + index)); /* mass in Msun - New correction: take into account the overdensity to get the halo mass*/
		  nhalos_cat++;
		}	    
	      }	
	    }
	  }
	  if(nhalos_cat>0) elem=fwrite(halo_cat,sizeof(Halo_t),nhalos_cat,fid_out);
	}  
	printf("Number of halos for M=%E: %ld\n",halo_mass, nhalos_R);fflush(0);	
      } /* matches else for when sigma is low */
      R=R/rhalo;
      halo_mass=(4.0/3.0)*PI*global_rho_m*(pow(R,3)+pow(R/rhalo,3))/2.0;  /* This is a new change - it seems this mass fits theory better... */
    } /* R cycle for excursion set */
    printf("--------------- Number of halos without subgrid: %ld -----------------------\n\n",nhalos_z);fflush(0);  

    
    /************************************************/
    /************************************************/
    /************************************************/
    if(global_use_sgrid==1){
      while (halo_mass >= global_halo_Mmin) {
	printf("\n\n----------------------------------------- subgrid halo mass=%E -----------------------------------\n", halo_mass); fflush(0);
	mass_aux=(4.0/3.0)*PI*global_rho_m*pow(R,3);
	nhalos_R=0;
	dm=(4.0/3.0)*PI*global_rho_m*(pow(R,3)-pow(R/rhalo,3));  /* mass interval */
	sigma_aux = sigma(R);   
	dndm_av = mass_function_ST(redshift, mass_aux)*dm*pow(global_dx_halo,3);  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume */
	hbias = Bias(redshift, sigma_aux);
	for (i = 0; i < global_N_halo; i++) {
	  for (j = 0; j < global_N_halo; j++) {
	    nhalos_cat=0;
	    for (p = 0; p < global_N_halo; p++) {
	      index = i * global_N_halo * global_N_halo + j * global_N_halo + p;                      
	      tempv=map_in_f[index]*growth*hbias;
	      Pnhalos = gsl_ran_poisson(rpoisson,dndm_av);
	      nn=(int)round(Pnhalos*(1.0+tempv));
	      for (ii = 0; ii < nn; ii++) {
		halo_cat[nhalos_cat].x=(int)i;
		halo_cat[nhalos_cat].y=(int)j;
		halo_cat[nhalos_cat].z=(int)p;
		halo_cat[nhalos_cat].Mass=halo_mass;
		nhalos_R++;
		nhalos_cat++;
		nhalos_z++;
	      }
	    }
	    if(nhalos_cat>0) elem=fwrite(halo_cat, sizeof(Halo_t), nhalos_cat, fid_out);
	  }
	}
	R = R/rhalo; 
	halo_mass=(4.0/3.0)*PI*global_rho_m*(pow(R,3)+pow(R/rhalo,3))/2.0;  /* This is a new change - it seems this mass fits theory better... */
	printf("Number of halos for M=%E: %ld\n",halo_mass, nhalos_R);fflush(0);	
      }        
      printf("----------------------------- Total number of halos including subgrid: %ld ---------------------- \n\n", nhalos_z);   fflush(0);
    } /* subgrid cycle */
    rewind(fid_out);
    elem=fwrite(&nhalos_z,sizeof(long int),1,fid_out);
    fclose(fid_out);       
        
    } /* ends box cycle */

    /************************************************/
  } /* ends redshift cycle */
      
  fftwf_destroy_plan(pc2r);
  fftwf_destroy_plan(pr2c);
  exit(0);    

}


