
/******************************************************************************************************************
SimFast21
Desciption: This routines determines collapsed halos positions and masses
from the linear density field. The barrier condition is given by Seth & Thormen modified version
of the critical overdensity.
Halos with mass smaller than cell size determined using a subgrid method
Main issue: excursion set assumes spherical halos - the approach breaks down when halo radius is close to cell size. Below Rmin we only use the density fluctuation in each cell to predict number of halos with mass lower than what is given by Rmin (which can be several cell sizes).
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

  double dndm_av, hbias, tempv, mass_dx;
  long int nhalos_z;  /* number of halos for a box z */
  long int nhalos_R; /* number of halos for radius R */
  double Pnhalos;
  int nhalos_cat;  /* number of halos to be written to catalog for a given slice */
  int nn;
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
  int flag, neg, zer;
  long int ncells1D;
  long int ind;
  long int index, indice;
  long int *halo_ind;
  char  *flag_halo;
  double deltaCMZ;
  double halo_mass;
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

  /* Memory allocations */
  if(!(halo_cat=(Halo_t *) fftwf_malloc(global_N_halo*global_N_halo*sizeof(Halo_t)))) {    /* get memory for float map */
    printf("Memory problem: halo_cat\n");
   exit(1);
  } 
  if(!(map_in_f=(float *) fftwf_malloc(global_N3_halo*sizeof(float)))) {    /* get memory for float map */
    printf("Memory problem - map_in_f...\n");
   exit(1);
  } 
  if(!(mass=(float *) fftwf_malloc(global_N3_halo*sizeof(float)))) {    /* get memory for float map */
    printf("Memory problem - mass...\n");
   exit(1);
  } 
  if(!(map_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
    printf("Memory problem - map_in...\n");
    exit(1);
  }

  /* FFTW plan */
  if(!(pr2c=fftwf_plan_dft_r2c_3d(global_N_halo, global_N_halo, global_N_halo, map_in_f , map_in, FFTW_ESTIMATE))) { 
    printf("Problem - fftwf_plan...\n");
    exit(1);
  }    

  /* Reading dark matter density box */
  printf("Reading dark matter density box\n");fflush(0);
  sprintf(dens_filename, "%s/delta/delta_z0_N%ld_L%.1f.dat", argv[1],global_N_halo, global_L/global_hubble);  
  fid_in=fopen(dens_filename,"rb");	/* second argument contains name of input file */
  if (fid_in==NULL) {
    printf("Error reading delta density file %s - check if the file exists...\n",dens_filename); 
    exit (1);
  }  
  elem=fread(map_in_f,sizeof(float),global_N3_halo,fid_in);    
  fclose(fid_in);
  fftwf_execute(pr2c); /* FFT of density box */
  
  if(!(map_in_aux = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
    printf("Memory problem - map_in_aux...\n");
    exit(1);
  }
  if(!(flag_halo=(char *) malloc(global_N3_halo*sizeof(char)))) {  
    printf("Memory problem - flag_halo...\n");
    exit(1);
  }
  if(!(halo_ind=(long int *) malloc(global_N3_halo/200*sizeof(long int)))) { 
    printf("Memory problem - halo_ind...\n");
    exit(1);
  }
  //  printf("smooth_f: %ld\n",smooth_factor);  
  //  printf("Preparing FFT from filtered k-space to real space...\n");fflush(0);
  if(!(pc2r=fftwf_plan_dft_c2r_3d(global_N_halo, global_N_halo, global_N_halo,(fftwf_complex *) map_in_aux,(float *)map_in_aux , FFTW_MEASURE ))) { 
    printf("Problem - FFTW...\n");
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


  mass_dx=global_rho_m*pow(global_dx_halo,3)*pow(10.0,global_halo_dlm); /* sets minimum halo mass for standard cycle (adds "a bit" to dx) */
  printf("Mininum halo mass for excursion set formalism: %E\n",mass_dx);fflush(0);
  if(mass_dx > global_halo_Mmin && global_use_sgrid==0)
    printf("Warning - minimum mass for box resolution is larger than halo_Mmin: you need to use subgrid to achieve that minimum mass.\n");

  
  
  /****************************************************/
  /****************************************************/
  /********  Redshift cycle ********/
  printf("Redshift cycle...\n");fflush(0);
  for(redshift=zmax;redshift>(zmin-dz/10);redshift-=dz){
    
    printf("\n+++++++++++++++++++++++++++ z = %f +++++++++++++++++++++++++++++++\n",redshift);fflush(0);
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
  
      growth=getGrowth(redshift);       /* get growth factor to evolve density field */
      sprintf(halo_filename, "%s/Halos/halo_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,(global_L/global_hubble));
      fid_out=fopen(halo_filename,"wb");
      if(fid_out==NULL){printf("\n Catalog file error\n"); exit (1);}
      nhalos_z=0;  /* counts total halos for one box */
      elem=fwrite(&nhalos_z,sizeof(long int),1,fid_out); /* reserve space for number of halos at beginning of file */    
      R=global_halo_Rmax;
      halo_mass=(4.0/3.0)*PI*global_rho_m*pow(R,3);  /* starting mass */
 

      /**************** start halo cycle for one box ****************/
      while(halo_mass>mass_dx && halo_mass>=global_halo_Mmin){   

	dndm_av=mass_function_ST(redshift, halo_mass)*halo_mass*(pow(10.0,global_halo_dlm)-1)*global_L3;  /* this is the expected average number of halos in the box within the target mass bin */
//	printf("dndm_av: %E\n",dndm_av);
	if(dndm_av >-log(1-0.1)) {  /* if there is more than a 10% (0.1) Poisson probability of finding a halo of this mass, go on... */
	  R2=R*R;   
	  sigma_aux=sigma(R);
	  deltaCMZ=growth*deltaFilter(sigma_aux,growth);

	  
	  /*************************************************/	
	  /*cycle to find halos for a given mass/radius start */
	  /* evolves density field, FFT and multiplies by window function */
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
		*((fftwf_complex *)map_in_aux + (i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p))= *((fftwf_complex *)map_in + (i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p))*W_filter(kk*R)*global_dx_halo*global_dx_halo*global_dx_halo/global_L3*growth;  /* This is the smoothed delta on a scale R.  Take care of all norms now */	      
	      }
	    }
	  }	
	  //	  printf("Executing FFT to real space...\n");fflush(0);
	  /* Executes FFT */
	  fftwf_execute(pc2r);  /* Inline transform - careful with output format... */
	  
	  
	  /* cycle start: go over box and find halos - checks for overlap */
	  /************************************************/
	  nhalos_R=0; /* halo count for one size R */	
	  ncells1D=(long int)(R/global_dx_halo); /* size of halo radius in cell width */
	  if(ncells1D==0)ncells1D=1;      
	  //	printf("Checking for halos and overlap...\n");fflush(0);
	  for(i=0;i<global_N_halo;i++){ 
	    nhalos_cat=0; /* counts number of halos written to catalogue - write catalogue to disk for each box slice instead of full box to avoid memory problems */
	    for(j=0;j<global_N_halo;j++){    
	      for(p=0;p<global_N_halo;p++){
		index=i*2*global_N_halo*(global_N_halo/2+1)+j*2*(global_N_halo/2+1)+p;  /* note the inline padding... */
		if(*((float *)map_in_aux + index)>deltaCMZ){  /* found a halo... */
		  flag=0;
		  ind=0;
		  /* check for halo overlap */
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
			    if(flag_halo[a*global_N_halo*global_N_halo+b*global_N_halo+c]!=1){ /* if find overlaps, breaks cycle */
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
		    for(indice=0;indice<ind;indice++)flag_halo[halo_ind[indice]]=1;  /* flag all mass cells that were used for this halo */
		    halo_cat[nhalos_cat].x=i;
		    halo_cat[nhalos_cat].y=j;
		    halo_cat[nhalos_cat].z=p;
		    halo_cat[nhalos_cat].Mass=halo_mass; /* mass in Msun */
		    //		  halo_cat[nhalos_cat].Mass=halo_mass*(1+*((float *)map_in_aux + index)); /* mass in Msun - New correction: take into account the overdensity to get the halo mass*/
		    nhalos_cat++;
		  }	    
		} /* cycle that finds halos */	
	      }
	    }
	    if(nhalos_cat>0) elem=fwrite(halo_cat,sizeof(Halo_t),nhalos_cat,fid_out);
	  } /* matches cycle over i */  
	  printf("Number of halos for M=%E: %ld\n",halo_mass, nhalos_R);fflush(0);
	  /* cycle end: go over box and find halos - checks for overlap */
	  /************************************************/
	/* cycle to detect halos for a given mass/radius terminates */
	/*************************************************/	

	} /* matches if to check if there is a low probability */
	halo_mass=halo_mass/pow(10.0,global_halo_dlm);  /* next mass value */
	R=pow(3/4.0/PI/global_rho_m*halo_mass,1.0/3);
      } /* mass cycle for excursion set ends */
      printf("\n--------------- Number of halos without subgrid: %ld -----------------------\n\n",nhalos_z);fflush(0);  

      
      /*******************************************/
      /* start cycle to check halos with cell size */
      halo_mass=global_rho_m*pow(global_dx_halo,3);
      nhalos_R=0;
      if(halo_mass >= global_halo_Mmin) {	
	R=pow(3/4.0/PI/global_rho_m*halo_mass,1.0/3);
	sigma_aux=sigma(R);
	deltaCMZ=growth*deltaFilter(sigma_aux,growth);
	for (i = 0; i < global_N_halo; i++) {
	  for (j = 0; j < global_N_halo; j++) {
	    nhalos_cat=0;
	    for (p = 0; p < global_N_halo; p++) {
	      index = i * global_N_halo * global_N_halo + j * global_N_halo + p;
	      if(flag_halo[index]==0) {
		if(map_in_f[index]*growth>deltaCMZ) {
		  flag_halo[index]=1;
		  halo_cat[nhalos_cat].x=i;
		  halo_cat[nhalos_cat].y=j;
		  halo_cat[nhalos_cat].z=p;
		  halo_cat[nhalos_cat].Mass=halo_mass;
		  nhalos_R++;
		  nhalos_cat++;
		  nhalos_z++;
		}
	      }
	    }
	    if(nhalos_cat>0) elem=fwrite(halo_cat, sizeof(Halo_t), nhalos_cat, fid_out);
	  }
	}
	printf("Number of halos found in cell size cycle: %ld\n",nhalos_R);
      }
      
      halo_mass=halo_mass/pow(10.0,global_halo_dlm);
      R=pow(3/4.0/PI/global_rho_m*halo_mass,1.0/3);

      /*** subgrid cycle for same z ***/
      /************************************************/
      if(global_use_sgrid==1){
	neg=0;
	zer=0;
	while (halo_mass >= global_halo_Mmin) {
	  //	printf("\n----------------------------------------- subgrid halo mass=%E -----------------------------------\n", halo_mass); fflush(0);
	  nhalos_R=0;
	  sigma_aux = sigma(R);   
	  dndm_av = mass_function_ST(redshift, halo_mass)*halo_mass*(pow(10.0,global_halo_dlm)-1)*pow(global_dx_halo,3);  /* mass_function_ST in 1/Msun/(Mpc/h)^3 - comoving volume */
	  hbias = Bias(redshift, sigma_aux);
	
	  for (i = 0; i < global_N_halo; i++) {
	    for (j = 0; j < global_N_halo; j++) {
	      nhalos_cat=0;
	      for (p = 0; p < global_N_halo; p++) {
		index = i * global_N_halo * global_N_halo + j * global_N_halo + p;
		if(flag_halo[index]==0)  {
		  tempv=map_in_f[index]*growth*hbias;
		  Pnhalos = gsl_ran_poisson(rpoisson,dndm_av); /* number of halos for a given mass given by Poisson distribution */
		  nn=(int)round(Pnhalos*(1.0+tempv));
		  if (tempv < -1) neg++;
		  for (ii = 0; ii < nn; ii++) {
		    halo_cat[nhalos_cat].x=i;
		    halo_cat[nhalos_cat].y=j;
		    halo_cat[nhalos_cat].z=p;
		    halo_cat[nhalos_cat].Mass=halo_mass;
		    nhalos_R++;
		    nhalos_cat++;
		    nhalos_z++;
		  }
		}
	      }
	      if(nhalos_cat>0) elem=fwrite(halo_cat, sizeof(Halo_t), nhalos_cat, fid_out);
	    }
	  }

	  halo_mass=halo_mass/pow(10.0,global_halo_dlm);
	  R=pow(3/4.0/PI/global_rho_m*halo_mass,1.0/3);

	// printf("Number of halos for M=%E: %ld\n",halo_mass, nhalos_R);fflush(0);	
	} /* ends subgrid cycle over halo masses */
	/*************************************/
	//	printf("Negative values: %d\n", neg);
	printf("----------------------------- Total number of halos including subgrid: %ld ---------------------- \n\n", nhalos_z);   fflush(0);
      } /* subgrid cycle ends */
      rewind(fid_out);
      elem=fwrite(&nhalos_z,sizeof(long int),1,fid_out);
      fclose(fid_out);       
        
    }
  } /* ends redshift cycle */
   /****************************************************/
  /****************************************************/
  /****************************************************/
     
  fftwf_destroy_plan(pc2r);
  fftwf_destroy_plan(pr2c);
  exit(0);    

}


