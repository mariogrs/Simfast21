
/*********************************************************************************************************
SimFast21
Description: adjust_halos - nonlinear correction of halo position using zeldovich approximation
*********************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>  

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){

  
  FILE *fid;
  long int i;
  float x,y,z;
  long int nhalos;
  long int indice;  
  float *map_veloc_realx;
  float *map_veloc_realy;
  float *map_veloc_realz;
  double redshift;
  Halo_t *halo_v;
  size_t elem;
  char fname[300];
  double zmin,zmax,dz;
  double growth;
 

  if(argc != 2) {
    printf("Nonlinear corrections to halo catalog using Zel'dovich approximation.\n");
    printf("usage: adjust_halos base_dir\n");
    printf("base_dir contains simfast21.ini\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  zmin=global_Zminsim;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;
 
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif

 
  /**********************************************************************************/
  /*************  Read velocity fields for nonlinear corrections  *************/

  printf("Reading Velocity...\n"); fflush(0);   
  if(!(map_veloc_realx=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in float */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_x_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,(global_L/global_hubble)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading X velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realx,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realy=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in float */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_y_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,(global_L/global_hubble)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Y velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realy,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realz=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in float */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_z_z0_N%ld_L%.1f.dat", argv[1],global_N_halo,(global_L/global_hubble)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Z velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realz,sizeof(float),global_N3_halo,fid);
  fclose(fid);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_veloc_realx,map_veloc_realy,map_veloc_realz,global_N_halo,global_N3_halo,global_L) private(i)
#endif  
  // linear displacement field
  for(i=0;i<(global_N3_halo);i++) {
    map_veloc_realx[i]*=(global_N_halo/global_L)*global_hubble;  /* correct for the fact that now the velocity in file is in Mpc not Mpc/h */
    map_veloc_realy[i]*=(global_N_halo/global_L)*global_hubble;
    map_veloc_realz[i]*=(global_N_halo/global_L)*global_hubble;
  }
  printf("Velocity Reading done...\n"); fflush(0);   

  
 /********************************************************************************************************/
 /******************************* Start redshift cycle over boxes ****************************************/

  printf("Redshift cycle...\n");fflush(0);
  for(redshift=zmax;redshift>(zmin-dz/10);redshift-=dz){
    
    printf("z = %f\n",redshift);fflush(0);
    growth=getGrowth(redshift)-getGrowth(300);
    
    /* read halo catalog */
    sprintf(fname, "%s/Halos/halo_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,(global_L/global_hubble));
    if((fid=fopen(fname,"rb"))==NULL){  
      printf("Halo file: %s does not exist... Check path or run get_halos for this configuration\n",fname);
      exit(1);
    } 
    elem=fread(&nhalos,sizeof(long int),1,fid);
    printf("Reading %ld halos...\n",nhalos);fflush(0);
    if(!(halo_v=(Halo_t *) malloc(nhalos*sizeof(Halo_t)))) { 
      printf("Problem - halo...\n");
      exit(1);
    }
    elem=fread(halo_v,sizeof(Halo_t),nhalos,fid);  
    fclose(fid);    
    
    //    printf("Adjusting halo positions...\n");fflush(0);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(nhalos,halo_v,map_veloc_realx,map_veloc_realy,map_veloc_realz,global_N_halo) private(i,indice,x,y,z)
#endif  
    for(i=0;i<nhalos;i++){
      x= halo_v[i].x;     
      y= halo_v[i].y;  
      z= halo_v[i].z;    
      indice=(long int)(((int)x)*global_N_halo*global_N_halo+((int)y)*global_N_halo+(int)z);
      x += (map_veloc_realx[indice]*growth);         
      y += (map_veloc_realy[indice]*growth);         
      z += (map_veloc_realz[indice]*growth);             
      halo_v[i].x=x;
      halo_v[i].y=y;
      halo_v[i].z=z;   
    }

    printf("Writing full non-linear halo catalog\n");fflush(0);
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.1f.dat.catalog",argv[1],redshift,global_N_halo,(global_L/global_hubble)); 
    if((fid=fopen(fname,"wb"))==NULL){  
      printf("\nError opening Halonl output catalog\n");
      return 0;
    }
    elem=fwrite(&nhalos,sizeof(long int),1,fid);
    elem=fwrite(halo_v,sizeof(Halo_t),nhalos,fid); 
    fclose(fid);
    
    free(halo_v);
    
  } /* ends redshift */
   
 exit(0);    
 
}
