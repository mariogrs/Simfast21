
/***************************************************************************************************
SimFast21
Name: get_densityfield										     
Description: This routines generates a Monte Carlo realization of the linear density field 
normalized by the matter power spectrum form the Eisenstein&Hu fitting formulae (or CAMB) at z=0
****************************************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>  
#include <fftw3.h>    
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
 
#include "Input_variables.h"
#include "auxiliary.h"

double *matrix_power[2];
long int iMax_ps; /* number of points actually read from the power spectrum file*/
double factor;


/* Calculates sigma8 directly from a file (for instance from CAMB output) */
double sig8_file(double *kv, double *Pv, int N) {
  int i;
  double dk,k,P,sig8;

  sig8=0;
  for(i=0;i<N-1;i++) {
    dk=kv[i+1]-kv[i];
    k=(kv[i]+kv[i+1])/2.0;
    P=(Pv[i]+Pv[i+1])/2;      
    sig8=sig8+k*k*P/(2*PI*PI)*W2(k*8.0)*dk;
  }

}


/*Power spectrum at z=0  */
double power_k(double k, tf_parms *tf){
  
  double  pkk;
 
  pkk=powerFunction(k,tf);  //power in (Mpc/h)^3
 
  return pkk;
}

long int search_kindex(double kk){

  long int i;

  for(i=0;i<iMax_ps-1;i++){
    if(matrix_power[0][i]<= kk && matrix_power[0][i+1]>=kk) break;
  }
  return i;

}

/*Function yielding the overdensity in Fourier space */
double deltaK(double k,int flag, tf_parms *tf){
  
  double variance;
  long int ind;
  
  if(k>0) {
    if(flag==1) {
      ind=search_kindex(k);
      variance = matrix_power[1][ind]+(((matrix_power[1][ind+1]-matrix_power[1][ind])/(matrix_power[0][ind+1]-matrix_power[0][ind]))*(k-matrix_power[0][ind]));
    } else variance = power_k(k,tf);
    variance=sqrt(variance/2.)*factor;
  } else variance=0.;
  
  return variance;
}




/* Main function, constructs the linear density field in real space */

int main(int argc, char **argv){

  long int i, j, p;                             
  long int indi,indj;
  double kk, aux_k, aux_pk;
  double sig8_old;
  float *map_f;
  size_t elem;
  fftwf_complex *map_in;
  fftwf_plan pc2r;
  FILE *file_power;
  FILE *file_out;
  char fname[256];
  DIR* dir;
  gsl_rng * r;
  tf_parms tf;
  int eSetCosm;
  
  if(argc != 2) {
    printf("Usage: get_densityfield work_dir\n");
    printf("work_dir - directory containing simfast21.ini\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif


  /*************************************************************************************/
  /*Alocation of memory for the vectors that will contain the density fields in real and k-space*/
  
  if(!(map_f=(float *) fftwf_malloc(global_N_halo*global_N_halo*global_N_halo*sizeof(float)))) {    /* get memory for float map */
    printf("Problem1...\n");
    exit(1);
  } 
 
  if(!(map_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)* global_N_halo*global_N_halo*(global_N_halo/2+1)))) {  
    printf("Problem2...\n");
    exit(1);
  }
  /* Prepare FFT to real space (use map_f as output) */
  /* the complex input map is just map_in (same structure) */
  if(!(pc2r=fftwf_plan_dft_c2r_3d(global_N_halo, global_N_halo, global_N_halo, map_in, map_f, FFTW_ESTIMATE))) { 
    printf("Problem...\n");
    exit(1);
  }   
 
  if(global_pk_flag==1){
  
    sprintf(fname, "%s/%s",argv[1],global_pk_filename); 
    if((file_power=fopen(fname,"r"))==NULL){  
      printf("\nThe CAMB Pk file cannot be open\n");
      exit(1);;
    }
    iMax_ps=0;
    while(!feof(file_power)){     
      elem=fscanf(file_power,"%lf  %lf\n",&aux_k, &aux_pk);  
      iMax_ps++;    
    }
    rewind(file_power);
    matrix_power[0]=(double *) malloc(iMax_ps*sizeof(double));
    matrix_power[1]=(double *) malloc(iMax_ps*sizeof(double));
    iMax_ps=0;
    while(!feof(file_power)){
      elem=fscanf(file_power,"%lf  %lf\n",&matrix_power[0][iMax_ps],&matrix_power[1][iMax_ps]);  
      iMax_ps++;                                              
    }
    if(matrix_power[0][iMax_ps-1] < global_dk*(global_N_halo/2.0)*sqrt(3.1)){
      printf("Error:k-range in matter power spectrum file not large enough\n"); 
      printf("%f  %f\n",matrix_power[0][iMax_ps-1],global_dk*(global_N_halo/2)*sqrt(3.1));
      exit(1); 
   }
    sig8_old=sig8_file(matrix_power[0],matrix_power[1],iMax_ps);
  }else{
    sig8_old=sig8(global_omega_m, global_omega_b, global_lambda);
  }  
  factor=global_sig8_new/sig8_old; 
  /* Construct density field in k-space*/


  /* gets random values (period must be larger than global_N_halo**3)
     parallelization requires to "carefully" initialize the RNG for each thread - leave it for now...
     maybe see SPRNG? */
  /* Select RNG */
  r = gsl_rng_alloc (gsl_rng_mt19937);
  //  r = gsl_rng_alloc (gsl_rng_ranlxd1);
  //  r = gsl_rng_alloc (gsl_rng_taus);
  //  r = gsl_rng_alloc (gsl_rng_gfsr4);
  gsl_rng_set (r, global_seed);
  for(i=0;i<global_N_halo*global_N_halo*(global_N_halo/2+1);i++) map_in[i]=gsl_ran_gaussian(r,1.0)+I*gsl_ran_gaussian(r,1.0);

  eSetCosm= Set_Cosmology(global_omega_m, global_omega_b, global_lambda, 0.0, &tf); /* sets cosmology for transfer function calculation at z=0 */

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N_halo,global_dk,map_in,global_pk_flag) private(i,indi,j,indj,p,kk) 
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
	map_in[i*global_N_halo*(global_N_halo/2+1)+j*(global_N_halo/2+1)+p]*=deltaK(kk, global_pk_flag, &tf);
      }
    }
  }
  gsl_rng_free (r);

  /*Introduce border symmetries */
  box_symmetriesf(map_in,global_N_halo);
    
  for(i=0;i<global_N_halo*global_N_halo*(global_N_halo/2+1);i++){
    map_in[i]=map_in[i]*sqrt(global_L3);
  }
 
  /* Execute FFT */
  fftwf_execute(pc2r);
  fftwf_destroy_plan(pc2r);
  free(map_in);     

  /* adds the missing factors to normalize */
  for(i=0;i<global_N_halo*global_N_halo*global_N_halo;i++){
    map_f[i]/=global_L3;
  }
 /* Create directory delta */
  sprintf(fname,"%s/delta",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating delta directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }
  sprintf(fname, "%s/delta/delta_z0_N%ld_L%d.dat", argv[1],global_N_halo, (int)(global_L));  
  file_out=fopen(fname,"wb");
  if (file_out==NULL) {
    printf("\n Problem opening the delta output file:%s\n",fname); 
    exit (1);
  }
  elem=fwrite(map_f,sizeof(float),global_N3_halo,file_out);  
  fftwf_free(map_f);    
  fclose(file_out);

  exit(0);    
  
}

