
/************************************************************************************************************
SimFast21
Module: get_HIIbubbles
Description: Generates a 3d box with global_N_smooth^3 cells witho 0s (neutral) and 1s (ionized)
Input: the base directory where we find the parameter file simfast21.ini
Output: creates base_dir/Ionization if needed
For further details see: 
M. G. Santos, L. Ferramacho, M. B. Silva, A. Amblard, A. Cooray, MNRAS 2010, http://arxiv.org/abs/0911.2219 
*************************************************************************************************************/

#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>   /* header for complex numbers in c */
#include <fftw3.h>     /* headers for FFTW library */
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Input_variables.h"
#include "auxiliary.h"

#define bfactor 1.1  /* value by which to divide bubble size R */
#define xHlim 0.999 /* cutoff limit for bubble calculation */
#define FFTWflag FFTW_MEASURE  /* PATIENT is too slow... */
//#define FFTWflag FFTW_PATIENT  /* PATIENT is too slow... */

int main(int argc, char *argv[]) {
  
  char fname[300];
  FILE *fid;
  DIR* dir;
  size_t elem;
  long int ii,ij,ik, ii_c, ij_c, ik_c, a, b, c;   
  long int ncells_1D;
  long int i,j,p,indi,indj,ind;
  int flag_bub,iz;
  double redshift,tmp;
  double kk;
  double neutral,*xHI;
  float *halo_map, *top_hat_r, *density_map,*bubblef, *bubble;  
  fftwf_complex *halo_map_c, *top_hat_c, *collapsed_mass_c, *density_map_c, *total_mass_c, *bubble_c;
  fftwf_plan pr2c1,pr2c2,pr2c3,pr2c4,pc2r1,pc2r2,pc2r3;
  double zmin,zmax,dz;
  double R;
  
  
  if(argc == 1 || argc > 5) {
    printf("Generates boxes with ionization fraction for a range of redshifts\n");
    printf("usage: get_HIIbubbles base_dir [zmin] [zmax] [dz]\n");
    printf("base_dir contains simfast21.ini and directory structure\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  if(argc > 2) {
    zmin=atof(argv[2]);
    if (zmin < global_Zminsim) zmin=global_Zminsim;
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
    zmin=global_Zminsim;
    zmax=global_Zmaxsim;
    dz=global_Dzsim;
  }
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);fflush(0);
#endif
  /* Create directory Ionization */
  sprintf(fname,"%s/Ionization",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Ionization directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }    
  sprintf(fname,"%s/Output_text_files",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Output_text_files directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }  
  
  
  
 /* Memory allocation - we could do some of the FFTs inline... */
 /*************************************************************/
  
  /* density_map mass */ 
  if(!(density_map=(float *) fftwf_malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem1...\n");
    exit(1);
  }
  if(!(density_map_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem2...\n");
    exit(1);
  }
  if(!(pr2c1=fftwf_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, density_map, density_map_c, FFTWflag))) { 
    printf("Problem3...\n");
    exit(1);
  }  
  /* halo_map mass */ 
  if(!(halo_map=(float *) fftwf_malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem4...\n");
    exit(1);
  }
  if(!(halo_map_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem5...\n");
    exit(1);
  }
  if(!(pr2c2=fftwf_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, halo_map, halo_map_c, FFTWflag))) { 
    printf("Problem6...\n");
    exit(1);
  }  
  /* total mass */
  if(!(total_mass_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem7...\n");
    exit(1);
  }
  if(!(pc2r1=fftwf_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, total_mass_c, density_map, FFTWflag))) { 
    printf("Problem8...\n");
    exit(1);
  }    
  /* collapsed mass */
  if(!(collapsed_mass_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem9...\n");
    exit(1);
  }
  if(!(pc2r2=fftwf_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, collapsed_mass_c, halo_map, FFTWflag))) { 
    printf("Problem10...\n");
    exit(1);
  }  
  /* top hat window */
  if(!(top_hat_r=(float *) fftwf_malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem11...\n");
    exit(1);
  }
  if(!(top_hat_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem12...\n");
    exit(1);
  }
  if(!(pr2c3=fftwf_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, top_hat_r, top_hat_c, FFTWflag))) { 
    printf("Problem13...\n");
    exit(1);
  } 
  /* bubble boxes */
  if(!(bubble=(float *) fftwf_malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem14...\n");
    exit(1);
  }
  if(!(bubble_c=(fftwf_complex *) fftwf_malloc(global_N_smooth*global_N_smooth*(global_N_smooth/2+1)*sizeof(fftwf_complex)))) {
    printf("Problem15...\n");
    exit(1);
  }
  if(!(bubblef=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem16...\n");
    exit(1);
  }
  if(!(pr2c4=fftwf_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, bubble, bubble_c, FFTWflag))) { 
    printf("Problem17...\n");
    exit(1);
  } 
  if(!(pc2r3=fftwf_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, bubble_c, bubble, FFTWflag))) { 
    printf("Problem18...\n");
    exit(1);
  }  
  if(!(xHI=(double *) malloc((int)((zmax-zmin)/dz+2)*sizeof(double)))) {
    printf("Problem19...\n");
    exit(1);
  }

  /****************************************************/
  /***************** Redshift cycle *******************/
  printf("Number of bubble sizes: %d\n",(int)((log(global_bubble_Rmax)-log(2.*global_dx_smooth))/log(bfactor)));
  printf("Redshift cycle...\n");fflush(0);
  iz=0;
  neutral=0.;
  for(redshift=zmin;redshift<(zmax+dz/10) && (neutral < xHlim);redshift+=dz){    
    printf("z = %f\n",redshift);fflush(0);

    sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
    fid=fopen(fname,"rb");
    if (fid==NULL) {printf("\nError reading deltanl file... Check path or if the file exists..."); exit (1);}
    elem=fread(density_map,sizeof(float),global_N3_smooth,fid);
    fclose(fid);
    
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_smooth, density_map,global_rho_m, global_dx_smooth,bubblef) private(i)
#endif
    for(i=0;i<(global_N3_smooth);i++){
      density_map[i]=(1.0+density_map[i])*global_rho_m*global_dx_smooth*global_dx_smooth*global_dx_smooth; /* total mass in 1 cell */
      bubblef[i]=0.0;
    }
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
    fid=fopen(fname,"rb");
    if (fid==NULL) {printf("\nError reading %s file... Check path or if the file exists...",fname); exit (1);}
    elem=fread(halo_map,sizeof(float),global_N3_smooth,fid);
    fclose(fid);

    /* Quick fill of single cells before going to bubble cycle */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_smooth,halo_map,density_map,global_eff,bubblef) private(i,tmp)
#endif
    for(i=0;i<global_N3_smooth;i++) {
      if(halo_map[i]>0.) {
	if(density_map[i]>0.) 
	  tmp=(double)halo_map[i]*global_eff/density_map[i];
	else tmp=1.0;
      }else tmp=0.;
      if(tmp>=1.0) bubblef[i]=1.0; else bubblef[i]=tmp;
    }
    /* FFT density and halos */
    fftwf_execute(pr2c1);    
    fftwf_execute(pr2c2);

    /************** going over the bubble sizes ****************/
    R=global_bubble_Rmax;    /* Maximum bubble size...*/
    while(R>=2*global_dx_smooth){ 
    
      printf("R= %lf\n", R);fflush(0);    
      //      printf("Filtering halo and density boxes...\n");fflush(0);
    
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(collapsed_mass_c,halo_map_c,total_mass_c,density_map_c,global_N_smooth,global_dk,R) private(i,j,p,indi,indj,kk)
#endif
      for(i=0;i<global_N_smooth;i++) {
	if(i>global_N_smooth/2) {
	  indi=-(global_N_smooth-i);
	}else indi=i;      
	for(j=0;j<global_N_smooth;j++) {
	  if(j>global_N_smooth/2) {
	    indj=-(global_N_smooth-j);
	  }else indj=j;
	  for(p=0;p<=global_N_smooth/2;p++) {
	    kk=global_dk*sqrt(indi*indi+indj*indj+p*p);
	    total_mass_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=density_map_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]*W_filter(kk*R);
	    collapsed_mass_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=halo_map_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]*W_filter(kk*R);
	  }
	}
      }
      
      fftwf_execute(pc2r1);     
      fftwf_execute(pc2r2); 
      
      flag_bub=0;
      
      //      printf("Starting to find and fill bubbles...\n");fflush(0);

      /* signal center of bubbles */      
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(halo_map,density_map,bubble,global_N_smooth,global_eff,flag_bub) private(ii,ij,ik,ind)
#endif  
      for(ii=0;ii<global_N_smooth;ii++){
	for(ij=0;ij<global_N_smooth;ij++){
	  for(ik=0;ik<global_N_smooth;ik++){
	    ind=ii*global_N_smooth*global_N_smooth+ij*global_N_smooth+ik;
	    if(halo_map[ind]>0.) {
	      if(density_map[ind]>0.) { 
		if((double)halo_map[ind]/density_map[ind]>=1.0/global_eff) {
		  flag_bub=1;
		  bubble[ind]=1.0;  	     	    	    
		}else bubble[ind]=0;
	      }else {
		flag_bub=1;
		bubble[ind]=1.0;  	     	    	    
	      }
	    }else bubble[ind]=0;
	  }
	}
      }
    
      /* generate spherical window in real space for a given R */
      if(flag_bub>0){
	printf("Found bubble...\n");fflush(0);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(top_hat_r,R,global_dx_smooth,global_N_smooth) private(i,j,p)
#endif  
	for(i=0;i<global_N_smooth;i++){
	  for(j=0;j<global_N_smooth;j++){
	    for(p=0;p<global_N_smooth;p++){ 
	      if(sqrt(i*i+j*j+p*p)*global_dx_smooth<=R || sqrt(i*i+(j-global_N_smooth)*(j-global_N_smooth)+p*p)*global_dx_smooth<=R  || sqrt(i*i+(j-global_N_smooth)*(j-global_N_smooth)+(p-global_N_smooth)*(p-global_N_smooth))*global_dx_smooth <=R || sqrt(i*i+(p-global_N_smooth)*(p-global_N_smooth)+j*j)*global_dx_smooth<=R || sqrt((i-global_N_smooth)*(i-global_N_smooth)+j*j+p*p)*global_dx_smooth<=R ||  sqrt((i-global_N_smooth)*(i-global_N_smooth)+(j-global_N_smooth)*(j-global_N_smooth)+p*p)*global_dx_smooth<=R ||  sqrt((i-global_N_smooth)*(i-global_N_smooth)+(j-global_N_smooth)*(j-global_N_smooth)+(p-global_N_smooth)*(p-global_N_smooth))*global_dx_smooth<=R ||  sqrt((i-global_N_smooth)*(i-global_N_smooth)+j*j+(p-global_N_smooth)*(p-global_N_smooth))*global_dx_smooth<=R ) {
		top_hat_r[i*global_N_smooth*global_N_smooth+j*global_N_smooth+p]=1.0;
	      }else top_hat_r[i*global_N_smooth*global_N_smooth+j*global_N_smooth+p]=0.0;	  
	    }
	  }
	}
	/* FFT bubble centers and window */ 
	fftwf_execute(pr2c3);
	fftwf_execute(pr2c4);
      
	/* Make convolution */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(bubble_c,top_hat_c,global_N_smooth) private(i)
#endif
	for(i=0;i<global_N_smooth*global_N_smooth*(global_N_smooth/2+1);i++) {
	  bubble_c[i]*=top_hat_c[i];
	}
	fftwf_execute(pc2r3);     
	
	/* after dividing by global_N3_smooth, values in bubble are between 0 (neutral)and global_N3_smooth */     
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(bubble,bubblef,global_N3_smooth) private(i)
#endif
	for (i=0; i<global_N3_smooth; i++){
	  bubble[i]/=global_N3_smooth;
	  if (bubble[i]>0.2) bubblef[i]=1.0; /* neutral should be zero */
	} 
	
      } /* ends filling out bubbles in box for R */
    
      R/=bfactor;  
    } /* ends R cycle */
 
    /* just to check smallest bubbles through older method */
    printf("Going to smaller R cycle...\n"); fflush(0);
    while(R>global_dx_smooth){ /* removed = global_dx_smooth since this was checked at the beginning */ 
  
      printf("R= %lf\n", R);fflush(0); 
      flag_bub=0;
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(collapsed_mass_c,halo_map_c,total_mass_c,density_map_c,global_N_smooth,global_dx_smooth,global_dk,R) private(i,j,p,indi,indj,kk)
#endif
      for(i=0;i<global_N_smooth;i++) {
	if(i>global_N_smooth/2) {
	  indi=-(global_N_smooth-i);
	}else indi=i;      
	for(j=0;j<global_N_smooth;j++) {
	  if(j>global_N_smooth/2) {
	    indj=-(global_N_smooth-j);
	  }else indj=j;
	  for(p=0;p<=global_N_smooth/2;p++) {
	    kk=global_dk*sqrt(indi*indi+indj*indj+p*p);
	    total_mass_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=density_map_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]*W_filter(kk*R);
	    collapsed_mass_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]=halo_map_c[i*global_N_smooth*(global_N_smooth/2+1)+j*(global_N_smooth/2+1)+p]*W_filter(kk*R);
	  }
	}
      }
      fftwf_execute(pc2r1); 
      fftwf_execute(pc2r2); 
   
      /* fill smaller bubbles in box */
      ncells_1D=(long int)(R/global_dx_smooth);    
      //      printf("Starting to find and fill bubbles...\n");fflush(0);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(halo_map,density_map,global_eff,global_N_smooth,global_dx_smooth,R,ncells_1D,bubblef,flag_bub) private(ii_c,ij_c,ik_c,ii,ij,ik,a,b,c,ind)
#endif
      for(ii_c=0;ii_c<global_N_smooth;ii_c++){
	for(ij_c=0;ij_c<global_N_smooth;ij_c++){
	  for(ik_c=0;ik_c<global_N_smooth;ik_c++){
	    ind=ii_c*global_N_smooth*global_N_smooth+ij_c*global_N_smooth+ik_c;
	    if(halo_map[ind]>0.) {
	      if(!(density_map[ind]>0.) || ((double)halo_map[ind]/density_map[ind]>=1.0/global_eff)) {
		flag_bub=1;
		for(ii=-(ncells_1D+1);ii<=ncells_1D+1;ii++){
		  a=check_borders(ii_c+ii,global_N_smooth);
		  for(ij=-(ncells_1D+1);ij<=ncells_1D+1;ij++){
		    if(sqrt(ii*ii+ij*ij)*global_dx_smooth <= R){
		      b=check_borders(ij_c+ij,global_N_smooth);
		      for(ik=-(ncells_1D+1);ik<=ncells_1D+1;ik++){
			c=check_borders(ik_c+ik,global_N_smooth);
			if(sqrt(ii*ii+ij*ij+ik*ik)*global_dx_smooth <= R){
			  bubblef[a*global_N_smooth*global_N_smooth+b*global_N_smooth+c]=1.0;  	     
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      if(flag_bub>0){printf("Found bubble...\n");fflush(0);}

      R/=bfactor;
    } /* ends small bubbles R cycle */

    neutral=0.;
    for (i=0; i<global_N3_smooth; i++){
      neutral+=1.0-bubblef[i];
    }
    neutral/=global_N3_smooth;
    printf("neutral fraction=%lf\n",neutral);fflush(0);
    xHI[iz]=neutral;
    sprintf(fname, "%s/Ionization/xHII_z%.3f_eff%.2lf_N%ld_L%.0f.dat",argv[1],redshift,global_eff,global_N_smooth,global_L); 
    if((fid = fopen(fname,"wb"))==NULL) {
      printf("Error opening file:%s\n",fname);
      exit(1);
    }
    elem=fwrite(bubblef,sizeof(float),global_N3_smooth,fid);  
    fclose(fid);
    iz++;
  } /* ends redshift cycle */


  /* z cycle for neutral>=xHlim */
  while(redshift<(zmax+dz/10)) {
    printf("z(>%f) = %f\n",xHlim,redshift);fflush(0);
    xHI[iz]=1.0;
    sprintf(fname, "%s/Ionization/xHII_z%.3f_eff%.2lf_N%ld_L%.0f.dat",argv[1],redshift,global_eff,global_N_smooth,global_L); 
    if((fid = fopen(fname,"wb"))==NULL) {
      printf("Error opening file:%s\n",fname);
      exit(1);
    }
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(bubblef,global_N3_smooth) private(i)
#endif
    for(i=0;i<global_N3_smooth;i++) bubblef[i]=0.0;
    elem=fwrite(bubblef,sizeof(float),global_N3_smooth,fid);  
    fclose(fid);
    iz++;
    redshift+=dz;
  }    

  sprintf(fname, "%s/Output_text_files/zsim.txt",argv[1]);
  if((fid = fopen(fname,"a"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  for(redshift=zmax;redshift>(zmin-dz/10);redshift-=dz) fprintf(fid,"%f\n",redshift); /* first line should be highest redshift */
  fclose(fid);
  sprintf(fname, "%s/Output_text_files/x_HI_eff%.2lf_N%ld_L%.0f.dat",argv[1],global_eff,global_N_smooth,global_L);
  if((fid = fopen(fname,"a"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  for(i=iz-1;i>=0;i--) fprintf(fid,"%lf\n",xHI[i]); /* first line should be highest redshift */
  fclose(fid);
  
  free(xHI);
  free(bubblef);
  fftwf_free(top_hat_r);
  fftwf_free(top_hat_c);
  fftwf_free(collapsed_mass_c);
  fftwf_free(density_map);
  fftwf_free(density_map_c);
  fftwf_free(halo_map);
  fftwf_free(halo_map_c);
  fftwf_free(total_mass_c);
  fftwf_free(bubble);
  fftwf_free(bubble_c);


  exit(0);
}
