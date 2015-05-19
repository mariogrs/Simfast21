
/************************************************************************************************************
SimFast21
convert_deltanl
Description: Converts density boxes (full dark matter simulations, e.g. non-linear) to Simfast21 format (the density_nl boxes).
M. G. Santos
*************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char *argv[]) {
  
  char fname[300];
  FILE *fid1, *fid2;
  size_t elem;
  int i,j,p, nx, ny, nz;
  double redshift;
  float *densitynl_in, *densitynl_out;
  double aver;
  
  if(argc!=4) {
    printf("\nConverts density filed (non-linear) from Mellema to Simfast21 format\n");
    printf("Assumes that the size of the input box is given by global_L and size cell in one dimension is given by global_N_smooth.\n");
    printf("usage: convert_deltanl.x  work_dir  input_density_file  redshift\n\n");
    exit(1);
  }  

  get_Simfast21_params(argv[1]);

  redshift=atof(argv[3]);
  sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
  fid1=fopen(fname,"rb");
  if(fid1!=NULL) {
    printf("File:%s already exists - exiting...\n",fname); 
    exit(1);
  }
  if((fid1=fopen(fname,"wb"))==NULL){  
    printf("\nError opening output nl_density file... Check path...\n");
    exit(1);
  } 

  /* read deltanl simulation box to convert */
  if((fid2=fopen(argv[2],"rb"))==NULL){  
    printf("Original density box file: %s does not exist...\n",argv[2]);
    exit(1);
  } 
  elem=fread(&nx,sizeof(int),1,fid2);
  elem=fread(&ny,sizeof(int),1,fid2);
  elem=fread(&nz,sizeof(int),1,fid2);
  printf("box cell number: (%d,%d,%d)\n",nx,ny,nz); 
  if((nx!=ny) || (nx!=nz) || (ny!=nz)) {printf("Error: box not cubic! Exiting...\n");exit(1);}
  if(nx != global_N_smooth) {
    printf("Error: box resolution different from global_N_smooth in the simfast21.ini! Exiting...\n");
    exit(1);
  }

  /* density nl */ 
  if(!(densitynl_in=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Mem problem...\n");
    exit(1);
  }
  if(!(densitynl_out=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Mem problem...\n");
    exit(1);
  }
  printf("Converting file...\n");
  elem=fread(densitynl_in,sizeof(float),global_N3_smooth,fid2);
  fclose(fid2);
  aver=0.;
  for(i=0;i<nx; i++) {
    for(j=0;j<nx; j++) {
      for(p=0;p<nx; p++) {
	aver+=(double)densitynl_in[p*nx*nx+j*nx+i];
      }
    }
  }
  for(i=0;i<nx; i++) {
    for(j=0;j<nx; j++) {
      for(p=0;p<nx; p++) {
	densitynl_out[i*nx*nx+j*nx+p]=(float)(((double)(densitynl_in[p*nx*nx+j*nx+i])-aver)/aver); /* swap indices */
      }
    }
  }
  elem=fwrite(densitynl_out,sizeof(float),global_N3_smooth,fid1);                    
  fclose(fid1);
    
  exit(0);
}


