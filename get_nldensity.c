
/**************************************************************************************************************
SimFast21
Output: delta(z) with nonlinear corrections
**************************************************************************************************************/

/* --------------Includes ----------------------------------------- */
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Input_variables.h"
#include "auxiliary.h"


int main(int argc, char **argv){
 
  FILE *fid;
  long int ind;
  long int i,j,p,x,y,z;
  size_t elem;
  double redshift;
  double growth;
  double growth_in;
  float *map_veloc_realx;
  float *map_veloc_realy;
  float *map_veloc_realz;
  float *map_in;
  float *map_out;
  float *map_out2;
  char fname[300];
  double zmin,zmax,dz;

  if(argc == 1 || argc > 5) {
    printf("Generates non-linear density boxes for a range of redshifts\n");
    printf("usage: get_nldensity base_dir [zmin] [zmax] [dz]\n");
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
  printf("Using %d threads\n",global_nthreads);
#endif

  printf("\nzmin=%E zmax=%E  dz=%E\n",zmin,zmax,dz); 
 
  /**********************************************************************************/
  /*************  Read velocity fields for nonlinear corrections  *************/

  printf("Reading velocity boxes...\n"); fflush(0);   
  if(!(map_veloc_realx=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem vx...\n");
    exit(1);
  }
  
  sprintf(fname, "%s/Velocity/vel_x_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading X velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realx,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realy=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem vy...\n");
    exit(1);
  }
  sprintf(fname, "%s/Velocity/vel_y_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Y velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realy,sizeof(float),global_N3_halo,fid);
  fclose(fid);
  if(!(map_veloc_realz=(float *) malloc(global_N3_halo*sizeof(float)))) {  /* get memory for output values in double */
    printf("Problem vz...\n");
    exit(1);
  }
  
  sprintf(fname, "%s/Velocity/vel_z_z0_N%ld_L%d.dat", argv[1],global_N_halo,(int)(global_L)); 
  fid=fopen(fname,"rb");	/* second argument contains name of input file */
  if (fid==NULL) {printf("\nError reading Z velocity file... Check path or if the file exists..."); exit (1);}
  elem=fread(map_veloc_realz,sizeof(float),global_N3_halo,fid);
  fclose(fid);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_veloc_realx,map_veloc_realy,map_veloc_realz,global_N_halo,global_N3_halo,global_L) private(ind)
#endif  
  //Calculo do campo de deslocamento linear
  for(ind=0;ind<global_N3_halo;ind++) {
    map_veloc_realx[ind]*=(global_N_halo/global_L);
    map_veloc_realy[ind]*=(global_N_halo/global_L);
    map_veloc_realz[ind]*=(global_N_halo/global_L);
  }
  printf("Velocity reading done...\n"); fflush(0);   

  if(!(map_out=(float *) malloc(global_N3_halo*sizeof(float)))) {  
    printf("Problem...\n");
    exit(1);
  }
  if(!(map_out2=(float *) malloc(global_N3_smooth*sizeof(float)))) {  
    printf("Problem...\n");
    exit(1);
  }
      
  if(!(map_in=(float *) malloc(global_N3_halo*sizeof(float)))) { 
    printf("Problem...\n");
    exit(1);
  }
  sprintf(fname, "%s/delta/delta_z0_N%ld_L%d.dat", argv[1],global_N_halo, (int)(global_L));  
  fid=fopen(fname,"rb");	
  if (fid==NULL){printf("\nError reading density file... Check if the file exists...\n"); exit (1);}
  elem=fread(map_in,sizeof(float),global_N3_halo,fid);
  fclose(fid);

  growth_in=getGrowth(300);
  
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(growth_in,global_N3_halo,map_in) private(ind)
#endif
  for(ind=0;ind<global_N3_halo;ind++) {
    map_in[ind]=map_in[ind]*growth_in+1.0;   
  }

  /***********************************************************************/
  /*********** Redshift cycle ********************************************/
  for(redshift=zmax;redshift>(zmin-dz/10.);redshift-=dz){
    
    sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
    if((fid=fopen(fname,"rb"))!=NULL){  
      printf("File %s already exists - skipping...\n",fname);
      fclose(fid);
    } else {
      growth=getGrowth(redshift)-growth_in;
      printf("Redshift: %lf\n",redshift);fflush(0);   
      
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_halo,map_out) private(ind)
#endif
      for(ind=0;ind<global_N3_halo;ind++) {
	map_out[ind]=-1.0;    
      }
      
    //#ifdef _OMPTHREAD_
    //#pragma omp parallel for shared(map_out,map_in,growth,global_N_halo) private(i,j,p,ind,x,y,z)
    //#endif 
      printf("Adjusting...\n");fflush(0);   
      for(i=0;i<global_N_halo;i++){
	for(j=0;j<global_N_halo;j++){
	  for(p=0;p<global_N_halo;p++){
	    ind=i*global_N_halo*global_N_halo+j*global_N_halo+p;
	    x=(long int)(i + map_veloc_realx[ind]*growth);
	    y=(long int)(j + map_veloc_realy[ind]*growth);
	    z=(long int)(p + map_veloc_realz[ind]*growth);   
	    x=check_borders(x,global_N_halo);
	    y=check_borders(y,global_N_halo);
	    z=check_borders(z,global_N_halo);
	    //#pragma omp critical
	    map_out[x*global_N_halo*global_N_halo+y*global_N_halo+z]+=map_in[ind];    
	  }
	}
      }
      printf("Smoothing...\n");fflush(0);   
      smooth_boxb(map_out, map_out2, global_N_halo, global_N_smooth);
      printf("Writing...\n");fflush(0);   
      sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_smooth,global_L); 
      if((fid=fopen(fname,"wb"))==NULL){  
	printf("\nError opening outpout nl_density file... Check path...\n");
	exit(1);
      } 
      elem=fwrite(map_out2,sizeof(float),global_N3_smooth,fid);                    
      fclose(fid);
    
      if(global_save_original_deltanl==1){
	sprintf(fname, "%s/delta/deltanl_z%.3f_N%ld_L%.0f.dat",argv[1],redshift,global_N_halo,global_L); 
	if((fid=fopen(fname,"wb"))==NULL){  
	  printf("\nError opening output nl_density file... Check path...\n");
	  exit(1);
	} 
	elem=fwrite(map_out,sizeof(float),global_N3_halo,fid);  
	fclose(fid);
      }
    }
      
  }

 
  exit(0);    

}
