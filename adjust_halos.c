
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
  long int i,j,p;
  long int x,y,z;
  long int nhalos;
  long int indice;  
  float *map_veloc_realx;
  float *map_veloc_realy;
  float *map_veloc_realz;
  float *halo_map;
  float *map_fcoll;
  double redshift;
  Halo_t *halo_v;
  size_t elem;
  char fname[300];
  double zmin,zmax,dz;
  double growth;
 

  if(argc==1 || argc > 5) {
    printf("Nonlinear corrections to halo catalog + generates halo mass box\n");
    printf("usage: adjust_halos base_dir [zmin] [zmax] [dz]\n");
    printf("base_dir contains simfast21.ini\n");
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
  printf("Using %d threads\n",global_nthreads);
#endif


 
  /**********************************************************************************/
  /*************  Read velocity fields for nonlinear corrections  *************/

  printf("Reading Velocity...\n"); fflush(0);   
  if(!(map_veloc_realx=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_x_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading X velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realx,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realy=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_y_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Y velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realy,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realz=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_z_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Z velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realz,sizeof(float),global_N3_halo,fid);
  fclose(fid);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_veloc_realx,map_veloc_realy,map_veloc_realz,global_N_halo,global_N3_halo,global_L) private(i)
#endif  
  //Calculo do campo de deslocamento linear
  for(i=0;i<(global_N3_halo);i++) {
    map_veloc_realx[i]*=(global_N_halo/global_L);
    map_veloc_realy[i]*=(global_N_halo/global_L);
    map_veloc_realz[i]*=(global_N_halo/global_L);
  }
  printf("Velocity Reading done...\n"); fflush(0);   
  
  if(!(halo_map=(float *) malloc(global_N3_smooth*sizeof(float)))) {
    printf("Problem...\n");
    exit(1);
  }
  

  printf("global_save_nl_halo_cat=%d\n",global_save_nl_halo_cat);
   
 /********************************************************************************************************/
 /******************************* Start redshift cycle over boxes ****************************************/

  printf("Redshift cycle...\n");fflush(0);
  for(redshift=zmax;redshift>(zmin-dz/10);redshift-=dz){
    
    printf("z = %f\n",redshift);fflush(0);
    growth=getGrowth(redshift)-getGrowth(300);
    
    /* read halo catalog */
    sprintf(fname, "%s/Halos/halo_z%.3f_N%ld_L%.0f.dat.catalog",argv[1],redshift,global_N_halo,global_L);
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
    
    
    printf("Adjusting...\n");fflush(0);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(nhalos,halo_v,map_veloc_realx,map_veloc_realy,map_veloc_realz,global_N_halo) private(i,indice,x,y,z)
#endif  
    for(i=0;i<nhalos;i++){
      x= halo_v[i].x;     
      y= halo_v[i].y;  
      z= halo_v[i].z;    
      indice=(long int)(x*global_N_halo*global_N_halo+y*global_N_halo+z);
      x += (long int)(map_veloc_realx[indice]*growth);         
      y += (long int)(map_veloc_realy[indice]*growth);         
      z += (long int)(map_veloc_realz[indice]*growth);             
      x=check_borders(x,global_N_halo);
      y=check_borders(y,global_N_halo);
      z=check_borders(z,global_N_halo);
      halo_v[i].x=x;
      halo_v[i].y=y;
      halo_v[i].z=z;   
    }
    
    printf("Collapsed box...\n");fflush(0);
    get_collapsed_mass_box(halo_map,halo_v, nhalos);
    if(global_save_nl_halo_cat==1){
      sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.0f.dat.catalog",argv[1],redshift,global_N_halo,global_L); 
      if((fid=fopen(fname,"wb"))==NULL){  
	printf("\nError opening Halonl output catalog\n");
	return 0;
      }
      elem=fwrite(halo_v,sizeof(Halo_t),nhalos,fid); 
      fclose(fid);
    }
    free(halo_v);
    
    if(global_use_fcoll==1){
      printf("Reading fcoll...\n");fflush(0);
      sprintf(fname, "%s/Halos/halo_fcoll_mass_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_halo,global_L);
      fid=fopen(fname,"rb");
      if(fid==NULL) printf("File:%s is missing - proceeding without fcoll mass...\n",fname);
      else {
	if(!(map_fcoll=(float *) malloc(global_N3_halo*sizeof(float)))) {  
	  printf("Mem Problem...\n");
	  exit(1);
	}
	elem=fread(map_fcoll,sizeof(float),global_N3_halo,fid);
	fclose(fid);
	printf("Adjusting  halos...\n");fflush(0);
	//#ifdef _OMPTHREAD_
	//#pragma omp parallel for shared(global_N_halo,map_poiss,map_veloc_realx,map_veloc_realy,map_veloc_realz,global_smooth_factor,halo_map,growth,global_N_smooth) private(i,j,p,indice,x,y,z)
	//#endif  
	for(i=0;i<global_N_halo;i++){
	  for(j=0;j<global_N_halo;j++){
	    for(p=0;p<global_N_halo;p++){
	      indice=i*global_N_halo*global_N_halo+j*global_N_halo+p;
	      if(map_fcoll[indice] > 100.) {  /* only do for nonzero masses */
		x = i+(long int)(map_veloc_realx[indice]*growth);         
		y = j+(long int)(map_veloc_realy[indice]*growth);         
		z = p+(long int)(map_veloc_realz[indice]*growth);             
		x=check_borders(x,global_N_halo);
		y=check_borders(y,global_N_halo);
		z=check_borders(z,global_N_halo);
		x=x/global_smooth_factor;
		y=y/global_smooth_factor;
		z=z/global_smooth_factor;
		//#pragma omp critical
		halo_map[x*global_N_smooth*global_N_smooth+y*global_N_smooth+z]+=map_fcoll[indice];
	      }
	    }
	  }
	}
	free(map_fcoll);
      }
    } /* ends fcoll*/

    printf("Writing box...\n");fflush(0);
    sprintf(fname, "%s/Halos/halonl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
    if((fid=fopen(fname,"wb"))==NULL){  
      printf("\nError opening Halonl output box\n");
      return 0;
    }
    elem=fwrite(halo_map,sizeof(float),global_N3_smooth,fid);
    fclose(fid);
    
  } /* ends redshift */
   
 exit(0);    
 
}
