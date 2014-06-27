
/* 
SimFast21
Generate epsilonX divided by n (see eq. 22 in paper Santos et al 2008)
output units in joules/s 
*/

#include <stdlib.h>
#include <stdio.h>
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <fftw3.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Input_variables.h"


int nzsfr;
double * zbox;
double * fHI;
double * fHeI;
double * fHeII;

const char * species_name[3]={"fHI","fHeI","fHeII"};
const double niovern[]={0.927,0.073,0.073};
const double Enth[]={13.6,24.59,54.42}; /* Ionization Energy threshold
					   for H,HeI,HeII */
/* The following values are the minimum and maximum values and step for the energy/frequency integration used in create_kernval_table */
double Enumax=5000.; /* eV */
double Enumin=60.; /* eV  - if it is less than Enth uses values above instead for the minimum in the integral */
int nEnu=1500;  /* number of Enu points for integration - don't like
		    this - do log integration maybe?  */


typedef struct {
  double xmin;
  double xmax;
  double dx;
  double *y;
} table;


typedef struct {
    double zmin;
    double zmax;
    double dz;
    double *r;
} r_table;

typedef struct {
    double rmin;
    double rmax;
    double dr;
    double *z;
} z_table;

double max(double val1, double val2);
double * produceKernel(double z, double zp, double zmax, int species); 
double emis(double z, double zp, double (* eb)(double));
double * convertfloat2doublearrandcombine(float * arr1, float * arr2);
double * convertfloat2doublearr(float * arr);
double verner_cross_section(double, int);
void convolve(fftw_complex * sfr_fft, double * kernel);
fftw_complex * MyFFT3D(double * sfr);


double get_int_r(double z); 
void create_rz_table(double zmin, double zmax, int n, r_table *rtab, z_table *ztab); 
void create_n_table(double zmin, double zmax, int n, table *tab);
void create_kernval_table(double zmin,double zmax, int n, table *tab);

double get_kernval(double z, int species);
double get_r(double z, r_table *rtab); 
double get_z(double r, z_table *ztab);
double drdz(double z);
double get_n(double,int);
void clean_rz_table(r_table *rtab, z_table *ztab);
void clean_n_table(table *tab);
void clean_kernval_table(table *tab);

/* Table to convert r to z and opposite*/


r_table rtab;
z_table ztab;
table tablen[3];
table tableker[3];


int main(int argc, char * argv[]) {

  double ztabmax,zrmax,ztocomputemax;
  //  double aver1,aver2;
  FILE * file;
  long i,box;
  double * sfr,  * Epsilon, *sfra;
  float * sfrt;
  float * fHII;
  double * kernel;
  fftw_complex * sfr_fft;
  char fname[256];
  int nztab=100000; /* we do need that high for at least get_n and get_kernval*/
  int species,ind,nzbox;
  double rc,drc,JXconst;
  double dzhom=0.001;  /* for the homogeneous integration */
  double zmin,zmax,dz,z;
  DIR* dir;
  int numboxes;
  //  int negval;
  float tmp;


  /* Check for correct number of parameters*/
  if(argc!=2) {
    printf("Calculates epsilonX/n for X ray temperature calculation\n");
    printf("usage: epsilonXon base_dir\n");
    printf("base_dir contains simfast21.ini\n");
    exit(1);
  }  
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) {printf("Lya and xray use set to false - no need to calculate xrays\n");exit(0);}
  zmin=global_Zminsfr;
  zmax=global_Zmaxsim;
  dz=global_Dzsim;
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif


  /* Create directory xrays */
  sprintf(fname,"%s/xrays",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating xray directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  } 

  nzsfr=(int)((zmax-zmin)/dz)+1;
  zbox = (double*) malloc(sizeof(double)*nzsfr);
  fHI = (double*) malloc(sizeof(double)*nzsfr);
  fHeI = (double*) malloc(sizeof(double)*nzsfr);
  fHeII = (double*) malloc(sizeof(double)*nzsfr);
  Epsilon = (double*) malloc(sizeof(double)*global_N3_smooth);
  sfra = (double *) malloc(sizeof(double)*(nzsfr));
  fHII = (float*) malloc(sizeof(float)*global_N3_smooth);                             

  sprintf(fname,"%s/Output_text_files/zsim.txt",argv[1]);
  if((file = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  for (i=0;i<nzsfr;i++) 
    if(fscanf(file,"%lf",&(zbox[nzsfr-1-i]))!=1) { 
      printf("Error reading zsim.txt!\n");
      exit(1);
    }
  fclose(file);    
  sprintf(fname, "%s/Output_text_files/x_HI_eff%.2lf_N%ld_L%.0f.dat",argv[1],global_eff,global_N_smooth,global_L);
  if((file = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  /* The average fHI and average SFR have the same redshifts - its up to the other program to make sure this happens */
  for (i=0;i<nzsfr;i++) {
    if(fscanf(file,"%lf",&(fHI[nzsfr-1-i]))!=1) {  /* starts again with highest redshift */
      printf("Error reading average HI file:%s\n",fname);
      exit(1);
    }
  }
  fclose(file);
  /* Assume HeI ionization is the same as HI and no HeII ionization */
  for (i=0;i<nzsfr;i++) {
    fHeI[i]=fHI[i];
    fHeII[i]=1.-fHeI[i];
  }
  sprintf(fname,"%s/Output_text_files/sfrd_av_N%ld_L%.0f.dat",argv[1],global_N_smooth,global_L);
  if((file = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  /* sfr in Msun/Mpc^3*h^3/yr (comoving volume and proper time) */
  /* Unit conversions are done in create_kernval_table */
  for (i=0;i<nzsfr;i++) {
    if(fscanf(file,"%f %lf",&tmp,&(sfra[nzsfr-1-i]))!=2) {  /* starts again with highest redshift */
      printf("Error reading %s file\n",fname);
      exit(1);
    }
  }
  fclose(file);

    
  /* I think that this should be ~ numax/numin*z (taking into account the homogeneous part) but that is quite high...
     we're limited to the available values of the SFR from the simulation (could use homogeneous calculation up 
     to higher redshifts but it seems to be negligible...) */
  ztocomputemax=zbox[nzsfr-1]+global_Dzsim-0.0000001;

  /**************************************************/
  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(nzbox=nzsfr-1;nzbox >= 0;nzbox--) {
    
    /* our box is at nzbox */
    printf("\n\nztocompute: %f\n",zbox[nzbox]);fflush(0);
    sprintf(fname,"%s/xrays/EpsilonXon_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
    if((file = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping this redshift...\n",fname);
      fclose(file);
    }else {
      /* create interpolation tables */
      create_rz_table(zbox[nzbox], zbox[nzsfr-1]+global_Dzsim, nztab, &rtab, &ztab);
      create_n_table(zbox[nzbox],zbox[nzsfr-1]+global_Dzsim,nztab,tablen);
      create_kernval_table(zbox[nzbox],zbox[nzsfr-1]+global_Dzsim,nztab,tableker);
   
      if(ztab.rmin+global_flux_Rmax>ztab.rmax) zrmax=ztocomputemax;
      else {
	zrmax = get_z(global_flux_Rmax+ztab.rmin,&ztab); 
	if(zrmax>ztocomputemax) zrmax=ztocomputemax;
      }
      zrmax=zrmax/1.001; /* reduce by 0.1% to avoid zrmax in the limit of boxes... */
      printf("zrmax: %f\n",zrmax);
      /*Compute how many box we need to convolve*/
      /* zbox points to beginning of redshft of the box */
      /* check if global_flux_Rmax is too small - in that case use only the average SFR */
      if(global_flux_Rmax < global_L/global_N_smooth*1.5) numboxes=0; 
      else {
	for (i=0;zbox[i+nzbox]+global_Dzsim < zrmax;i++);
	numboxes=i+1;
      }
      /* maximum numboxes is nzsfr */
      
      /*Initialize some variables*/  
      memset(Epsilon,0,global_N3_smooth*sizeof(double));
      sprintf(fname,"%s/Ionization/xHII_z%.3f_eff%.2lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_eff,global_N_smooth,global_L);
      printf("Read ioni file %s\n",fname);fflush(0);
      if((file = fopen(fname,"r"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fread(fHII,sizeof(float),global_N3_smooth,file);  
      fclose(file);
        
      /*************************************/
      /***** inhomogeneous computation *****/
      /*************************************/
      printf("Inhomogeneous calculation from %f to %f - using %d boxes\n",zbox[nzbox],zrmax,numboxes);fflush(0);
      for(box=0;box<numboxes;box++) {
	printf("Read SFR %li at z=%5.3lf\n",box,zbox[box+nzbox]);fflush(0);
	/* SFR in comoving (h/Mpc)**3 */
	sprintf(fname,"%s/SFR/sfrd_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[box+nzbox],global_N_smooth,global_L);
	if((file = fopen(fname,"r"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	sfrt = (float*) malloc(sizeof(float)*global_N3_smooth);
	fread(sfrt,sizeof(float),global_N3_smooth,file);
	fclose(file);
	sfr = convertfloat2doublearr(sfrt);
	free(sfrt);
	sfr_fft=MyFFT3D(sfr);
	
	/* Compute the kernel*/
	if (zbox[box+nzbox]+global_Dzsim >= zrmax) ztabmax=zrmax; else ztabmax=zbox[box+nzbox]+global_Dzsim;
	for (species=0;species<3;species++) {
	  kernel = produceKernel(zbox[nzbox],zbox[nzbox+box],ztabmax,species);            
	  /* Convolve kernel with sfr  */
	  convolve(sfr_fft,kernel);            
	  /*not species specific, we assume ionization fraction is the same*/
	  /* We need the HII boxes for all redshifts requested */
	  /* NOTE: this is the ionization box not fHI!! */
	  /* Assume ionization HI and HeI are the same, assume no HeIII so fHeII=1-fHeI=1-fHI*/
	  /* convert ionization to fHI! */
	  /* Addition to Epsilon */
	  if (species < 2) {        
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Epsilon,kernel) private(i) 
#endif
	    for (i=0;i<global_N3_smooth;i++) Epsilon[i]+=kernel[i]*(1.-fHII[i])*niovern[species];
	}else {
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Epsilon,kernel) private(i) 
#endif
	    for (i=0;i<global_N3_smooth;i++) Epsilon[i]+=kernel[i]*fHII[i]*niovern[species];
	  }
	  free(kernel);
	}
	free(sfr_fft);
      }

     
      /***********************************/
      /***** homogeneous computation *****/
      /***********************************/
      printf("Homogeneous calculation from %f to %f...\n",zrmax,ztocomputemax);fflush(0);
      JXconst=0.; 
      for (ind=nzbox;zbox[ind]+global_Dzsim <= zrmax;ind++);    
      rc=get_r(zrmax,&rtab)-ztab.rmin;
      for (z=zrmax; (z<ztocomputemax) && (z<zbox[nzsfr-1]+global_Dzsim); z+=dzhom) {
	drc=drdz(z)*dzhom;
	if(z>zbox[ind]+global_Dzsim) ind++;
	JXconst+=drc*rc*rc*sfra[ind]*(get_kernval(z,0)*fHI[ind]*niovern[0]+get_kernval(z,1)*fHeI[ind]*niovern[1]+get_kernval(z,2)*fHeII[ind]*niovern[2]);
	rc=rc+drc;  /* comoving Mpc/h */
      }
      JXconst*=4.*PI/(global_L/global_N_smooth)/(global_L/global_N_smooth)/(global_L/global_N_smooth);  /* cancels the dr^3 in create_kernval_table */
      
      /*
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Epsilon,JXconst) private(i) 
#endif 
      for (i=0;i<global_N3_smooth;i++) Epsilon[i]+=JXconst;
      */
      printf("Value of the Constant part %le\n",JXconst);
      
      for (i=0;i<global_N3_smooth;i++) Epsilon[i]+=JXconst;
      
      /* Write out the results */
      /*
      aver1=0.;
      aver2=0.;
      negval=0;
      for(i=0;i<global_N3_smooth;i++) {
	aver1+=Epsilon[i];
	if(Epsilon[i]<0.) {negval++;}
	//	if(Epsilon[i]<0.) {Epsilon[i]=0.; negval++;}
	aver2+=Epsilon[i];
      }
      if(aver1<0.) {
	printf("Error: epsilon average is negative! %E\n",aver1/global_N3_smooth);
	exit(1);
      }
      printf("%d negative values in a total of %ld cells (%f %%)\n",negval,global_N3_smooth,1.*negval/global_N3_smooth*100.);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Epsilon,aver1,aver2) private(i) 
#endif
      for(i=0;i<global_N3_smooth;i++) {
	Epsilon[i]*=aver1/aver2;
      }
      */
     
      sprintf(fname,"%s/xrays/EpsilonXon_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
      if((file = fopen(fname,"wb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fwrite(Epsilon,sizeof(double),global_N3_smooth,file);
      fclose(file);
      
      clean_rz_table(&rtab, &ztab);
      clean_n_table(tablen);
      clean_kernval_table(tableker);
    }
  } /* ends redshift cycle */
    
  
  exit(0);
}





double * produceKernel(double z, double zp, double zmax, int species) {

  double * Kernel;
  double zpp;
  long int ix,iy,iz,ind1,ind2,izmax,izmin,rmin,rmax;

  Kernel = (double*) malloc(sizeof(double)*global_N3_smooth);
  memset(Kernel,0,sizeof(double)*global_N3_smooth);

  /* Note: normalization is done in create_kernval_table */
  rmin = (long int) (get_r(zp, &rtab)-ztab.rmin)*global_N_smooth/global_L;  /* units are Mpc/h */
  rmax = (long int) (get_r(zmax, &rtab)-ztab.rmin)*global_N_smooth/global_L;
 
  if (rmax <= rmin) {
    printf("Problem in the Kernel production, rmax <= rmin,\n N:%li z:%lf zp:%lf zmax:%lf rmin:%li rmax:%li\n",global_N_smooth,z,zp,zmax,rmin,rmax);
  }else { 
  for (ix=-1*rmax;ix<=rmax;ix++) 
    //    printf("ix: %ld\n",ix);fflush(0);
    for (iy=-1*rmax;iy<=rmax;iy++) if ((ix*ix+iy*iy) <= rmax*rmax) {
        izmax = (long int)round(sqrt(rmax*rmax-ix*ix-iy*iy));
        if ((ix*ix+iy*iy) > rmin*rmin) izmin = 0 ;
        else izmin=(long int)round(sqrt(rmin*rmin-ix*ix-iy*iy));
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Kernel,izmin,izmax,rmax,ztab,ix,iy) private(iz,ind1,ind2,zpp)
#endif
	for (iz=izmin;iz<izmax;iz++) 
	  if (((ix*ix+iy*iy+iz*iz) < rmax*rmax) &&
	      (ix*ix+iy*iy+iz*iz !=0)){
	    
	    zpp = get_z(sqrt((double)(ix*ix+iy*iy+iz*iz))*global_L/global_N_smooth+ztab.rmin,&ztab);
	    
	    ind1=((long)(ix < 0)*global_N_smooth+ix)+((long)(iy < 0)*global_N_smooth+iy)*global_N_smooth+iz*global_N_smooth*global_N_smooth;
	    ind2=((long)(ix < 0)*global_N_smooth+ix)+((long)(iy < 0)*global_N_smooth+iy)*global_N_smooth+(global_N_smooth-iz)*global_N_smooth*global_N_smooth;
	    
	    /* Now kernel is tabulated and interpolated since
	       its value is spherically symmetric, it's a good speed-up*/
	    Kernel[ind1]=get_kernval(zpp,species);
	    if (iz !=0) Kernel[ind2]=Kernel[ind1];
	  }else if ((ix*ix+iy*iy+iz*iz ==0)) {
	    zpp = get_z(0.5*global_L/global_N_smooth+ztab.rmin,&ztab);  /* just to fill the zero point */
	    Kernel[0]=get_kernval(zpp,species);
	  }
      }
  }

  return Kernel;
}



fftw_complex * MyFFT3D(double * sfr) {

  fftw_complex * out1;
  fftw_plan plan;

 
  out1 = (fftw_complex*) malloc(sizeof(fftw_complex)*global_N_smooth*global_N_smooth*(global_N_smooth/2+1));

#ifdef _OMPTHREAD_
  fftw_init_threads();
  fftw_plan_with_nthreads(global_nthreads);
#endif
  
  plan=fftw_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, sfr, out1, FFTW_ESTIMATE);
  fftw_execute(plan);
  free(sfr);

#ifdef _OMPTHREAD_
  fftw_cleanup_threads();
#endif

  return out1;
}



void convolve(fftw_complex * sfr_fft, double * kernel) {

  fftw_complex *out2;
  double tempr,tempi;
  long np,i;
  fftw_plan plan;

 

#ifdef _OMPTHREAD_
  fftw_init_threads();
  fftw_plan_with_nthreads(global_nthreads);
#endif
  
  out2 = (fftw_complex*) malloc(sizeof(fftw_complex)*global_N_smooth*global_N_smooth*(global_N_smooth/2+1));
  plan=fftw_plan_dft_r2c_3d(global_N_smooth, global_N_smooth, global_N_smooth, kernel, out2, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan); 


  /* Maybe I have to do something extra for the normalization
     of the Fourier transform with N^3

  + don't forget the geometrical factor*/

  np = global_N_smooth/2 +1 ;
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(out2,sfr_fft,np) private(i,tempr,tempi)
#endif
  for (i=0;i<global_N_smooth*global_N_smooth*np;i++) {
    tempr = sfr_fft[i][0]*out2[i][0] - sfr_fft[i][1]*out2[i][1];
    tempi = sfr_fft[i][1]*out2[i][0] + sfr_fft[i][0]*out2[i][1];
    out2[i][0] = tempr/global_N3_smooth;
    out2[i][1] = tempi/global_N3_smooth;
  }


  plan=fftw_plan_dft_c2r_3d(global_N_smooth, global_N_smooth, global_N_smooth, out2, kernel, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan); 
  free(out2);  

#ifdef _OMPTHREAD_
  fftw_cleanup_threads();
#endif

}



double * convertfloat2doublearrandcombine(float * arr1, float * arr2) {
  double * out;
  long int i;
  
 
  out = (double*) malloc(sizeof(double)*global_N3_smooth);


#ifdef _OMPTHREAD_
#pragma omp parallel for shared(out,arr1,arr2,global_N3_smooth) private(i)
#endif 
  for (i=0;i<global_N3_smooth;i++) out[i] = (double) arr1[i] + (double) arr2[i];

  return out;
}

double * convertfloat2doublearr(float * arr) {
  double * out;
  long int i;
  

  out = (double*) malloc(sizeof(double)*global_N3_smooth);


#ifdef _OMPTHREAD_
#pragma omp parallel for shared(out,arr,global_N3_smooth) private(i)
#endif 
  for (i=0;i<global_N3_smooth;i++) out[i] = (double) arr[i];

  return out;
}



double verner_cross_section(double Enu, int species) {

  /* I only put data for HI, HeI and HeII since they are the only
     relevant here 
     I multiply by 1e-22 to convert the Mb into m^2
  */

/*------------------------------------------------------------------------------*/
/*   Bytes Format  Units   Label        Explanations                            */ 
/*------------------------------------------------------------------------------*/
/*   1-  2  I2     ---     Z            Atomic number                         0 */  
/*   4-  5  I2     ---     N            Number of electrons                   1 */
/*   7- 15  E9.3   eV      E_th         Subshell ionization threshold energy  2 */   
/*  17- 25  E9.3   eV      E_max        Maximum energy                        3 */ 
/*  27- 35  E9.3   eV      E_0          Fit parameter                         4 */
/*  37- 45  E9.3   Mb      \sigma_0     Fit parameter                         5 */ 
/*  47- 55  E9.3   ---     y_a          Fit parameter                         6 */
/*  57- 65  E9.3   ---     P            Fit parameter                         7 */
/*  67- 75  E9.3   ---     y_w          Fit parameter                         8 */
/*  77- 85  E9.3   ---     y_0          Fit parameter                         9 */  
/*  87- 95  E9.3   ---     y_1          Fit parameter                        10 */
/*------------------------------------------------------------------------------*/


  const double partable[3][11]={{1,1,13.6, 50000,0.4298,5.475e4,32.88,2.963,0,0,0},
			     {2,2,24.59,50000,13.61,949.2,1.469,3.188,2.039,0.4434,2.136},
			     {2,1,54.42,50000,1.72,13690,32.88,2.963,0,0,0}};
  
  
  double E0,y0,y1,yw,P,Ya,sig0,Y,x,F;
  
  E0=partable[species][4];
  y0=partable[species][9];
  y1=partable[species][10];
  yw=partable[species][8];
  P=partable[species][7];
  Ya=partable[species][6];
  sig0=partable[species][5];

  x=Enu/E0-y0;
  Y=sqrt(x*x+y1*y1);
  F=((x-1)*(x-1)+yw*yw)*pow(Y,0.5*P-5.5)*pow((1+sqrt(Y/Ya)),-1.*P);

  
  return sig0*F*1e-22;

}



double get_n(double z, int species) {
  
  int ix;

  if((z<tablen[species].xmin) || (z>tablen[species].xmax)) {
    printf("!!ERROR1!!\n");
    printf("valeurs : %lf %lf %lf\n",z,tablen[species].xmin,tablen[species].xmax);
    exit(1);
  }

  ix=(int)round((z-tablen[species].xmin)/tablen[species].dx);

  return tablen[species].y[ix];
}



double get_kernval(double z, int species) {

  int ix;
  double zp;

  if((z<tableker[species].xmin) || (z>tableker[species].xmax)) {
    printf("!!ERROR1!!\n");
    printf("valeurs : %lf %lf %lf\n",z,tableker[species].xmin,tableker[species].xmax);
    exit(1);
  }

  zp=z-tableker[species].xmin;
  ix=(int)floor((zp)/tableker[species].dx);
  
  return (tableker[species].y[ix]*((ix+1)*tableker[species].dx-zp)+tableker[species].y[ix+1]*(zp-ix*tableker[species].dx))/tableker[species].dx;
}



void create_kernval_table(double zmin, double zmax, int nz, table *tab) {

  int  i, j, k;
  double *dl, *tau, *zv;
  double  z, tmp,dz,norm,normb,r,Enu;
  double dEnu[3];
 
  norm = (1+zmin)*(1+zmin)*pow(global_Enu0*(1+zmin),global_alphas+1)*
    global_L0/global_Enu0*global_L/global_N_smooth/global_hubble/Mpc2m/Mpc2m/4./PI*global_L*global_L/global_N_smooth/global_N_smooth; /* there is a dr^3 here but just one 1/h since the others cancel with r^2 below */  
  /* this converts JX into j/m^2(proper)/s/eV^2 (there is an 1/r^2 below) */

  dz= (zmax-zmin)/(nz-1.);
  for(j=0;j<3;j++) {
    tab[j].y=(double *) malloc(nz*sizeof(double));
    memset(tab[j].y,0,nz*sizeof(double));
  }
  for (j=0;j<3;j++) {
    tab[j].xmin=zmin;
    tab[j].xmax=zmax;
    tab[j].dx=dz;
  }
  zv=(double *) malloc(nz*sizeof(double));
  dl=(double *) malloc(nz*sizeof(double));
  
  /* calculate dl in tau integral to save time */
  for(i=0;i<nz;i++) {
    z=zmin+i*dz;
    zv[i]=z;
    dl[i]=drdz(z)*dz*Mpc2m/global_hubble/(1+z);   /* I assume that both nHI and
					 sigma_HI, etc are in meters!!
					 I add the division by h
					 ALmost forgot that dl need to
					 be in physical distance*/
  }

  tau=(double *) malloc(nz*sizeof(double));
  tau[0]=0.;  /* when z'=zmin */
  for (i=0;i<3;i++) dEnu[i]=(Enumax-max(Enth[i],Enumin))/(nEnu-1.);  /* if Enumin is larger than Enth use it instead */

  for (j=0;j<3;j++) {
    //    printf("kerntable - species: %d\n",j);fflush(0);
    for(k=0;k<nEnu;k++) {
      Enu = max(Enth[j],Enumin)+dEnu[j]*k;
      /* Calculate tau */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(tau,nz,Enu,zv,dl,zmin) private(i) 
#endif
      for(i=1;i<nz;i++) 
	  tau[i]=dl[i]*
	      (get_n(zv[i],0)*verner_cross_section((1+zv[i])/(1+zmin)*Enu,0)+
	       get_n(zv[i],1)*verner_cross_section((1+zv[i])/(1+zmin)*Enu,1)+
	       get_n(zv[i],2)*verner_cross_section((1+zv[i])/(1+zmin)*Enu,2));

      for (i=1;i<nz;i++) tau[i]+=tau[i-1];
      /* get kernel vector for interpolation */            
      tmp=(Enu-Enth[j])*verner_cross_section(Enu,j)*dEnu[j];  /* eV^2 * m^2 * JX units (above) */
      /*put part of \hat{epsilon}_nu in here to save time */

#ifdef _OMPTHREAD_
#pragma omp parallel for shared(tab,nz,j,tmp,tau,global_alphas,Enu,zv) private(i) 
#endif
      for(i=0;i<nz;i++) 
	tab[j].y[i]+=tmp*exp(-1.*tau[i]-(global_alphas+1)*log(Enu*(1.+zv[i])));
      /* prefer to make nz Multiplication than a pow*/			     
    }
  }
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(tab,norm,nz,zmin,dz,ztab,rtab) private(i,r,z,normb) 
#endif
  for (i=1;i<nz;i++) {  /* start at 1 to avoid r=0 */
    z=zmin+dz*i;
    r=get_r(z,&rtab)-ztab.rmin;   /* r in Mpc/h - cancels global_L/global_N_smooth*global_L/global_N_smooth above in norm */
    normb=norm/r/r*global_hubble*global_hubble*global_hubble;  /* I removed the (1+z)^-3 since SFRD is already in comoving (h/Mpc)**3 - this conversion should actually be done when writing the sfr files (makes more sense than here), also added the h^3 to remove h dependence */
    /* no simulation 216 dependent normalization value for SFR is used */
    tab[0].y[i]*=normb;
    tab[1].y[i]*=normb;
    tab[2].y[i]*=normb;
  }

  free(dl);
  free(tau);
  free(zv);
  /* Should this be zero?... */
  tab[0].y[0]=0;
  tab[1].y[0]=0;
  tab[2].y[0]=0;
 
}


void clean_kernval_table(table *tab) {
  free(tab[0].y);
  free(tab[1].y);
  free(tab[2].y);
}



void create_n_table(double zmin, double zmax, int n, table *tab) {

  int i,j,k;
  double z,dz;
  double nH0,nHe0;

  if(zmin < zbox[0]) {
    printf("Error create_n_table: zmin less than zbox[0]\n");fflush(0);
    exit(1);
  }
  if(zmax > zbox[nzsfr-1]+global_Dzsim) {
    printf("Error create_n_table: zmax too large: %f\n",zmax);fflush(0);
    exit(1);
  }
  for (i=0;i<3;i++)
    tab[i].y = (double*) malloc(sizeof(double)*n);

  dz = (zmax-zmin)/(n-1.);

  for (i=0;i<3;i++) {
    tab[i].xmin = zmin;
    tab[i].xmax = zmax;
    tab[i].dx = dz;
  }

  nH0=global_omega_b*global_hubble*global_hubble*3*H0*H0/8./PI/G*(1-YHe)/mH;
  nHe0=global_omega_b*global_hubble*global_hubble*3*H0*H0/8./PI/G*YHe/mHe;

  k=0;
  for (j=0;j<n;j++) {
    z=zmin+dz*j;
    if(z>=zbox[nzsfr-1]) {
      tab[0].y[j]=nH0*(1+z)*(1+z)*(1+z)*fHI[nzsfr-1];
      tab[1].y[j]=nHe0*(1+z)*(1+z)*(1+z)*fHeI[nzsfr-1];
      tab[2].y[j]=nHe0*(1+z)*(1+z)*(1+z)*fHeII[nzsfr-1];
    }else {      
      while(z>=zbox[k]+global_Dzsim) k++;
      if (k > nzsfr-2) {
	printf("Error in create_n_table - k: %d\n",k);fflush(0);
	exit(1);
      }
      tab[0].y[j]=nH0*(1+z)*(1+z)*(1+z)*(fHI[k+1]*(z-zbox[k])+fHI[k]*(zbox[k+1]-z))
	/(zbox[k+1]-zbox[k]);
      tab[1].y[j]=nHe0*(1+z)*(1+z)*(1+z)*(fHeI[k+1]*(z-zbox[k])+fHeI[k]*(zbox[k+1]-z))
	/(zbox[k+1]-zbox[k]);
      tab[2].y[j]=nHe0*(1+z)*(1+z)*(1+z)*(fHeII[k+1]*(z-zbox[k])+fHeII[k]*(zbox[k+1]-z))
	/(zbox[k+1]-zbox[k]);
    }
  }
  for(i=0;i<3;i++) {
    for(j=0;j<n;j++) {
      if(isnan(tab[i].y[j])) { 
	printf("Error in n!! %d %d\n\n",i,j);
	exit(1);
      }
    }
  }

}


void clean_n_table(table *tab) {
  free(tab[0].y);
  free(tab[1].y);
  free(tab[2].y);
}



/* get r(z) */
/* Comoving Mpc/h */
double get_int_r(double z) {
    
    double dz=0.001,r;
    int n,i;

    n=(int)(z/dz)+1;
    dz=z/n;

    r=0.;
    for(i=0; i<n;i++) {
        r+=drdz(i*dz+dz/2.);
    }
 
    return r*dz;
}

/* create table for r(z) and z(r) */
void create_rz_table(double zmin, double zmax, int n, r_table *rtab, z_table *ztab) {

    double dz,dr,z,r, rmin, rmax;
    int i,iz;

    (*rtab).r=(double *) malloc(n*sizeof(double));
    dz=(zmax-zmin)/(n-1);

    for(i=0; i<n; i++) {
        z=zmin+i*dz;
        (*rtab).r[i]=get_int_r(z);
    }
    (*rtab).dz=dz;
    (*rtab).zmin=zmin;
    (*rtab).zmax=zmax;
    rmin=(*rtab).r[0];
    rmax=(*rtab).r[n-1];
    (*ztab).rmin=rmin;
    (*ztab).rmax=rmax;

    (*ztab).z=(double *) malloc(n*sizeof(double));
    dr=(rmax-rmin)/(n-1);

    (*ztab).z[0]=zmin;
    (*ztab).z[n-1]=zmax;

    for(i=1; i<n-1; i++) {
        r=rmin+i*dr;
        for (iz=0;(*rtab).r[iz]<r;iz++);
        if(((*rtab).r[iz]-r) < (r-(*rtab).r[iz-1])) {
            (*ztab).z[i]=zmin+iz*dz;
        }else {
            (*ztab).z[i]=zmin+(iz-1)*dz;
        }
    }
    (*ztab).dr=dr;

}

void clean_rz_table(r_table *rtab, z_table *ztab) {
  free((*rtab).r);
  free((*ztab).z);
}



/* comoving Mpc/h  */
double get_r(double z, r_table *rtab) {

    int nz;

    if((z<(*rtab).zmin) || (z>(*rtab).zmax)) {
        printf("!!ERROR1!!\n");
        printf("valeurs : %lf %lf %lf\n",z,(*rtab).zmin,(*rtab).zmax);
        exit(1);
    }
    nz=(int)round((z-(*rtab).zmin)/(*rtab).dz);
    
    return (*rtab).r[nz];
    
}

double get_z(double r, z_table *ztab) {

    int nr;

    if((r<(*ztab).rmin) || (r>(*ztab).rmax)) {
        printf("!!ERROR2!!\n");
        printf("values : %lf %lf %lf\n",r,(*ztab).rmin,(*ztab).rmax);
        exit(1);
    }
    nr=(int)round((r-(*ztab).rmin)/(*ztab).dr);

    return (*ztab).z[nr];

}

/* dr/dz*dz is in comoving Mpc/h */

double drdz(double z) {

    return 2997.9/sqrt(global_omega_m*(1.+z)*(1.+z)*(1.+z)+global_lambda);  /* value in Mpc/h */
}



double max(double val1, double val2) {

  double mv;

  if(val1>=val2) mv=val1; else mv=val2;

  return mv;

}



