
/*************************************************************
SimFast21
Calculates xalpha (Lya - coupling see Santos et al 2008/2010)
*************************************************************/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Input_variables.h"


double nualpha; /* Hz */
double nubeta; /* Hz */
const int nmax=16;
const double frec[]={0.,0.,1.,0.,0.2609,0.3078,0.3259,0.3353,0.3410,0.3448,
                     0.3476,0.3496,0.3512,0.3524,0.3535,0.3543};
const double Jc0=5.552e-8;      // m^-2 s^-1 Hz^-1 sr^-1 at z=0   =27*A10*Tcmb0/16/pi/sigma_a/Tstar
int nzsfr;
#define nuLL 3.28932e15  /* Lya limit frequency */


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

void convolve(double * sfr, double * kernel, long N);
double * produceKernel(long int N, double z, double zp, double zmax, double (* eb)(double)); 
double emis(double z, double zp, double (* eb)(double));
double * convertfloat2doublearr(float * arr, long nele);

double ebG(double nu);
double get_int_r(double z); 
void create_rz_table(double zmin, double zmax, int n, r_table *rtab, z_table *ztab); 
double get_r(double z, r_table *rtab); 
double get_z(double r, z_table *ztab);
double drdz(double z);

/* Table to convert r to z and opposite*/


r_table rtab;
z_table ztab;


int main(int argc, char * argv[]) {


  DIR* dir;
  double zmax2,ztabmax,zrmax;
  double aver1,Jc;
  FILE * file,*fidtxt;
  double * zbox; 
  int box,numboxes,ind,nzbox;
  //  int negval;
  long i;
  double * sfr, * Jalpha;
  float * sfrt;
  double * kernel;
  char fname[256];
  int nztab=10000;
  double *sfra;
  double jalphaconst,z;
  double dz = 0.001;
  float tmp;
  //  double aver2;

  /* Check for correct number of parameters*/
  if (argc != 2) {
    printf("Usage : xalpha base_dir\n");
    exit(0);
  }
  
  get_Simfast21_params(argv[1]);
  if(global_use_Lya_xrays==0) {printf("Lya and xray use set to false - no need to calculate Lya\n");exit(0);}
  
#ifdef _OMPTHREAD_
  omp_set_num_threads(global_nthreads);
  printf("Using %d threads\n",global_nthreads);
#endif
 

  /* Create directory Lya */
  sprintf(fname,"%s/Lya",argv[1]);
  if((dir=opendir(fname))==NULL) {  
    printf("Creating Lya directory\n");
    if(mkdir(fname,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))!=0) {
      printf("Error creating directory!\n");
      exit(1);
    }
  }

  nzsfr=(int)((global_Zmaxsim-global_Zminsfr)/global_Dzsim)+1;
  zbox = (double*) malloc(sizeof(double)*nzsfr);
  sfra = (double *) malloc(sizeof(double)*(nzsfr));
  /*Initialize some variables*/  
  nualpha=nuLL*(1.-(1./4.));  /* Hz */
  nubeta=nuLL*(1.-(1./9.));  /* Hz */  
  Jalpha = (double*) malloc(sizeof(double)*global_N3_smooth);
  sfrt = (float*) malloc(sizeof(float)*global_N3_smooth);
  
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
  sprintf(fname,"%s/Output_text_files/sfrd_av_N%ld_L%.0f.dat",argv[1],global_N_smooth,global_L);
  if((file = fopen(fname,"r"))==NULL) {
    printf("Error opening file:%s\n",fname);
    exit(1);
  }
  /* sfr in Msun/Mpc^3*h^3/yr (comoving volume and proper time) */
  for (i=0;i<nzsfr;i++) {
    if(fscanf(file,"%f %lf",&tmp,&(sfra[nzsfr-1-i]))!=2) {  /* starts again with highest redshift */
      printf("Error reading %s file\n",fname);
      exit(1);
    }
    /* convert to baryon number and comoving density in 1/Mpc^3 */
    sfra[nzsfr-1-i]*=Msun/mbar*global_hubble*global_hubble*global_hubble; /* starts with highest redshift! */
  }
  fclose(file);
  sprintf(fname,"%s/Output_text_files/xa_av_N%ld_L%.0f.dat",argv[1],global_N_smooth,global_L); 
  if((fidtxt=fopen(fname,"a"))==NULL){
    printf("\nError opening output %s file...\n",fname); 
    exit(1);
  }  

  /**************************************************/
  /************** redshift cycle ********************/
  /**************************************************/
  for(nzbox=nzsfr-1;nzbox >= 0;nzbox--) {
    /* our box is at nzbox */
   
    printf("\n\nztocompute: %f\n",zbox[nzbox]);fflush(0);
    sprintf(fname,"%s/Lya/xalpha_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
    if((file = fopen(fname,"rb"))!=NULL) {
      printf("File:%s already exists - skipping this redshift...\n",fname);
      fclose(file);
    }else {
      zmax2 = (1.+zbox[nzbox])*32./27.-1.; /* (1-1/9)/(1-1/4)*/
      printf("Lya full max redshift: %f\n",zmax2);
      /**** The condition below is a bit strict - maybe we can use a lower zmax2 and still get the integral to converge?? ****/
      if(zmax2>=zbox[nzsfr-1]+global_Dzsim) {
	printf("Need higher redshift boxes for proper computation - zmax2: %f  zf: %f - setting zmax2 to simulation zf...\n",zmax2,zbox[nzsfr-1]+global_Dzsim);
	zmax2=zbox[nzsfr-1]+global_Dzsim-0.0000001;
      }
      create_rz_table(zbox[nzbox], zbox[nzsfr-1]+global_Dzsim , nztab, &rtab, &ztab);
      if(ztab.rmin+global_flux_Rmax>ztab.rmax) zrmax=zmax2;
      else {
	zrmax = get_z(global_flux_Rmax+ztab.rmin,&ztab); 
	if(zrmax>zmax2) zrmax=zmax2;
      }
      zrmax=zrmax/1.001; /* reduce by 0.1% to avoid zrmax in the limit of boxes... */
      printf("zrmax: %f\n",zrmax);
      /*Compute how many box we need to zrmax*/
      /* Note: zrmax for 50 Mpc/h is always less than zmax2 */
      /* !!! This also assumes that zrmax < zbox[nzbox-1] !!! */ 
      /* Note: zbox[] is the redshift for the box, at the beginning of the bin of width global_Dzsim */
      /* check if global_flux_Rmax is too small - in that case use only the average SFR */
      if(global_flux_Rmax < global_L/global_N_smooth*1.5) numboxes=0; else {
	for (i=0;zbox[nzbox+i]+global_Dzsim < zrmax;i++);
	numboxes=i+1;
      }
      /* maximum numboxes is nzsfr */
      //    printf("Test: %f %f %f   %f %f  %f %f\n",ztocompute,zmax2,zrmax,global_flux_Rmax,ztab.rmin,zbox[numboxes-1],global_Dzsim);  
      for(i=0;i<global_N3_smooth;i++) Jalpha[i]=0.;
      
      /*************************************/
      /***** inhomogeneous computation *****/
      /*************************************/
      /*Loop over the box to compute the contribution to J_alpha*/  
      printf("Using %d boxes\n",numboxes);fflush(0);
      for(box=0;box<numboxes;box++) {
	
	printf("Convolving box %d\n",box);fflush(0);
	/*Load the box in memory */
	sprintf(fname,"%s/SFR/sfrd_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox+box],global_N_smooth,global_L);
	if((file = fopen(fname,"r"))==NULL) {
	  printf("Error opening file:%s\n",fname);
	  exit(1);
	}
	fread(sfrt,sizeof(float),global_N3_smooth,file);
	fclose(file);
	sfr = convertfloat2doublearr(sfrt, global_N3_smooth);
	
	/* Compute the kernel */
	if (zbox[nzbox+box]+global_Dzsim >= zrmax) ztabmax=zrmax; else ztabmax=zbox[nzbox+box]+global_Dzsim;
	kernel = produceKernel(global_N_smooth,zbox[nzbox],zbox[nzbox+box],ztabmax,&ebG);
	
	/* Convolve kernel with sfr */
	convolve(sfr,kernel,global_N_smooth);
	
	/* Addition to Jalpha */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Jalpha,kernel) private(i) 
#endif
	for (i=0;i<global_N3_smooth;i++) Jalpha[i]+=kernel[i];
      
	free(kernel);     
      }
  
      /***********************************/
      /***** homogeneous computation *****/
      /***********************************/
      jalphaconst=0.; 
      for (ind=nzbox;zbox[ind]+global_Dzsim <= zrmax;ind++);    
      for (z=zrmax; (z<zmax2) && (z<zbox[nzsfr-1]+global_Dzsim); z+=dz) {
	if(z>zbox[ind]+global_Dzsim) ind++;
	jalphaconst+=drdz(z)/global_hubble*sfra[ind]*emis(zbox[nzbox],z,&ebG);  /* 1/h to convert drdz from Mpc/h to Mpc */
      }  
      jalphaconst*=(1.+zbox[nzbox])*(1.+zbox[nzbox])/4./PI*dz/Mpc2m/Mpc2m/y2s;  /* units in number/m^2/s/Hz/sr */
      
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Jalpha,jalphaconst) private(i) 
#endif 
      for (i=0;i<global_N3_smooth;i++) Jalpha[i]+=jalphaconst;
      
      Jc=Jc0*(1.+zbox[nzbox]);

      printf("Value of the Constant part %le\n",jalphaconst/Jc);

      /* Write out the results */    
      /*
      aver1=0.;
      aver2=0.;
      negval=0;
      for(i=0;i<global_N3_smooth;i++) {
	Jalpha[i]/=Jc;
	aver1+=Jalpha[i];
	if(Jalpha[i]<0.) {negval++;}
	//	if(Jalpha[i]<0.) {Jalpha[i]=0.; negval++;}
	aver2+=Jalpha[i];
      }
      printf("%d negative values in %ld cells (%f %%)\n",negval,global_N3_smooth,1.*negval/global_N3_smooth*100.);fflush(0);
      if(aver1<0.) {
	printf("Error: xalpha average is negative! %E\n",aver1/global_N3_smooth);
	exit(1);
      }
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Jalpha,aver1,aver2) private(i) 
#endif
      for(i=0;i<global_N3_smooth;i++) {
	Jalpha[i]*=aver1/aver2;
      }
      */
      aver1=0.;
      for(i=0;i<global_N3_smooth;i++) {
	Jalpha[i]/=Jc;
	aver1+=Jalpha[i];
      }
      for(i=0;i<global_N3_smooth;i++) sfrt[i]=(float)Jalpha[i];    
      sprintf(fname,"%s/Lya/xalpha_z%.3lf_N%ld_L%.0f.dat",argv[1],zbox[nzbox],global_N_smooth,global_L);
      if((file = fopen(fname,"wb"))==NULL) {
	printf("Error opening file:%s\n",fname);
	exit(1);
      }
      fwrite(sfrt,sizeof(float),global_N3_smooth,file);
      fclose(file);
      fprintf(fidtxt,"%f %.8E\n",zbox[nzbox],aver1/global_N3_smooth);
    }
    
  } /* ends redshift cycle */
  
  fclose(fidtxt);
  exit(0);

}




    /*
double cH_inv(double z);
const double A10=2.85e-15;      // s^-1
const double Tstar=0.068;       // K 
const double Tcmb0=2.725;       // K
const double sigma_a=1.105e-6 ; // m^2 Hz
const double h=0.7;
const double omm=0.27;
const double oml=0.73;
const double omb=0.044;
const double nuLL=3.28932e15;   // Hz
double nualpha; // Hz

  nH_bar=0.76*0.0223*1.12321e-5*pow(1.+z,3);  // cm^-3  - use n_p/n_n=6.5
  Tcmb=Tcmb0*(1.+z);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(NT,z,nualpha,fHI,nH_bar,rho_mat,temp,jalpha) private(i,tauGP,Salpha)
#endif
  for(i=0; i<NT; i++) {
    tauGP=sigma_a*cH_inv(z)/nualpha*(double)fHI[i]*nH_bar/0.01/0.01/0.01*(double)rho_mat[i]; // no units
    //    tauGP=sigma_a*cH_inv(z)/nualpha*(double)fHI[i]*nH_bar/0.01/0.01/0.01*1.0; // no units
    //    tauGP=sigma_a*cH_inv(z)/nualpha*1.0*nH_bar/0.01/0.01/0.01*1.0; // no units 
    Salpha=exp(-0.803*pow(temp,-2./3)*pow(1.e-6*tauGP,1./3));   // good for TK > 1K
    jalpha[i]=(float)Salpha*jalpha[i]/Jc;

    // cH^-1 in m 

double cH_inv(double z) {

  return 2997.9/h/sqrt(omm*pow(1.+z,3)+oml)*3.08568e22;

}
    */






double * produceKernel(long int N, double z, double zp, double zmax,
		       double (* eb)(double)) {

  double * Kernel;
  double zpp,dr,norm;
  long int ix,iy,iz,izmax,izmin,ind1,ind2;
  long int rmin,rmax;
 
  Kernel = (double*) malloc(sizeof(double)*N*N*N);
  memset(Kernel,0,sizeof(double)*N*N*N);
 
  /*  rmin = get_dis(z,zp) ;
  rmax = get_dis(z,zmax) ;*/
  rmin = (long int) (get_r(zp, &rtab)-ztab.rmin)*N/global_L;
  rmax = (long int) (get_r(zmax, &rtab)-ztab.rmin)*N/global_L;
  

  dr = global_L/N/global_hubble;   /* have dr^3 (from FFTW) over dr^2 (from area) - units of Mpc (because we divided by h) */
  norm = (1+z)*(1+z)/4./PI*dr/4./PI*Msun/Mpc2m/Mpc2m/y2s/mbar*global_hubble*global_hubble*global_hubble; /* also includes h^3 - correction factor for sfr which is in h^3/Mpc^3 */

  if (rmax <= rmin) {
    printf("Problem in the Kernel production, rmax <= rmin,\n N:%li z:%lf zp:%lf zmax:%lf rmin:%li rmax:%li\n",N,z,zp,zmax,rmin,rmax);
  }else {
  for (ix=-1*rmax;ix<=rmax;ix++)
    for (iy=-1*rmax;iy<=rmax;iy++) if ((ix*ix+iy*iy) <= rmax*rmax) {
	izmax = (long int)round(sqrt(rmax*rmax-ix*ix-iy*iy));
	if ((ix*ix+iy*iy) > rmin*rmin) izmin = 0 ;
	else izmin=(long int)round(sqrt(rmin*rmin-ix*ix-iy*iy));
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(Kernel,izmin,izmax,rmax,N,ztab,ix,iy) private(iz,ind1,ind2,zpp)
#endif
	for (iz=izmin;iz<izmax;iz++) 
	  if (((ix*ix+iy*iy+iz*iz) < rmax*rmax) &&
	      (ix*ix+iy*iy+iz*iz !=0)){
	  zpp = get_z(sqrt(ix*ix+iy*iy+iz*iz)*global_L/N+ztab.rmin,&ztab);
	  /*	  printf("%i %i %i\n",ix,iy,iz);*/
	  ind1=((long)(ix < 0)*N+ix)+((long)(iy < 0)*N+iy)*N+iz*N*N;
	  ind2=((long)(ix < 0)*N+ix)+((long)(iy < 0)*N+iy)*N+(N-iz)*N*N;
	  /*printf("%li %li %lf %lf\n",ind1,ind2,z,zpp);*/
	  
	  Kernel[ind1]=emis(z,zpp,eb)/((double)(ix*ix+iy*iy+iz*iz))*norm;
	  if (iz !=0) Kernel[ind2]=Kernel[ind1];
	  } else  if ((ix*ix+iy*iy+iz*iz ==0)) {
	    zpp = get_z(0.5*global_L/N+ztab.rmin,&ztab);
	    Kernel[0]=emis(z,zpp,eb)/((double)(0.5*0.5))*norm;
	  }
      }
  }

  return Kernel;
}


double emis(double z, double zp, double (*eb)(double)) {

  int n;
  double sumemis=0;
  double zrat;
  
  zrat = (1+zp)/(1+z);

  /* since zmax is not the same for all n
     I chose to cut nmax instead of zmax 
     but that should be the same
     therefore zmax should be zmax(2)*/
  for (n=2;(zrat<=(1-1./(n+1.)/(n+1.))/(1.-1./n/n)) && (n<nmax);n++)    
    sumemis+=eb(nuLL*(1-1./n/n)*zrat)*frec[n];

  return sumemis;
}


void convolve(double * sfr, double * kernel, long N) {

  fftw_complex * out1, *out2;
  double tempr,tempi;
  long np,i;
  fftw_plan plan;
 
  out1 = (fftw_complex*) malloc(sizeof(fftw_complex)*N*N*(N/2+1));

#ifdef _OMPTHREAD_
  fftw_init_threads();
  fftw_plan_with_nthreads(global_nthreads);
#endif
  
  /*printf("cal 1st plan\n");*/
  plan=fftw_plan_dft_r2c_3d(N, N, N, sfr, out1, FFTW_ESTIMATE);
  /*printf("cal 1st fft\n");*/
  fftw_execute(plan);
  /*printf("free sfr\n");*/
  free(sfr);

  out2 = (fftw_complex*) malloc(sizeof(fftw_complex)*N*N*(N/2+1));
  /*printf("cal 2nd plan\n");*/
  plan=fftw_plan_dft_r2c_3d(N, N, N, kernel, out2, FFTW_ESTIMATE);
  /*printf("cal 2nd fft\n");*/
  fftw_execute(plan);
  /*printf("destroy plan\n");*/
  fftw_destroy_plan(plan); 


  /* Maybe I have to do something extra for the normalization
     of the Fourier transform with N^3

  + don't forget the geometrical factor*/

  np = N/2 +1 ;
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(out2,out1,N,np) private(i,tempr,tempi)
#endif
  for (i=0;i<N*N*np;i++) {
    tempr = out1[i][0]*out2[i][0] - out1[i][1]*out2[i][1];
    tempi = out1[i][1]*out2[i][0] + out1[i][0]*out2[i][1];
    out2[i][0] = tempr/N/N/N;
    out2[i][1] = tempi/N/N/N;
  }

  free(out1);

  /*printf("cal 3rd plan\n");*/
  plan=fftw_plan_dft_c2r_3d(N, N, N, out2, kernel, FFTW_ESTIMATE);
  /*printf("cal 3rd fft\n");*/
  fftw_execute(plan);
  /*printf("destroy plan\n");*/
  fftw_destroy_plan(plan); 
  free(out2); 

#ifdef _OMPTHREAD_
  fftw_cleanup_threads();
#endif

}


double * convertfloat2doublearr(float * arr, long nele) {
  double * out;
  long i;

  out = (double*) malloc(sizeof(double)*nele);
  for (i=0;i<nele;i++) out[i] = (double) arr[i];

  return out;
}




/* generic function for number photons/baryon/freq.(Hz) */
/* Assumes A/nu^0.9 SED with 20000 Lyman photons/baryon */
/* We can change A and index 0.9 (A is degenerate with SFR efficiency */
double ebG(double nu) {

  return global_A_Lya*pow(nu,-global_alpha_Lya);

}

/* get r(z) in Mpc/h*/
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
        (*rtab).r[i]=get_int_r(z);  /* Mpc/h */
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


/* in Mpc/h */
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


/* r in Mpc/h */
double get_z(double r, z_table *ztab) {

    int nr;

    if((r<(*ztab).rmin) || (r>(*ztab).rmax)) {
        printf("!!ERROR2!!\n");
        printf("valeurs : %lf %lf %lf\n",r,(*ztab).rmin,(*ztab).rmax);
        exit(1);
    }
    nr=(int)round((r-(*ztab).rmin)/(*ztab).dr);

    return (*ztab).z[nr];

}

/* dr/dz in comoving Mpc/h */

double drdz(double z) {

    return 2997.9/sqrt(global_omega_m*(1.+z)*(1.+z)*(1.+z)+global_lambda);  /* value in Mpc/h */


}

