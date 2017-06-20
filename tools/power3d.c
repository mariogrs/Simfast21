
/*
 Power spectrum calculator
Mario G. Santos (2014)
*/

#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.1415926535897932385E0


#ifdef _OMPTHREAD_
int nthreads;
#endif

#define SQR(x) ((x)*(x))
double wind(int j,int k,int p,int Nx,int Ny,int Nz, int flag);



/*********************************************************************/
int main(int argc, char **argv){

  long int N1, N2, N3;       /* number of pixels on the side */
  int bsize;
  double dl, L1, L2, L3, k;     /* Mpc/h */
  long int i,j,p,indi, indj, nk, *nps, indk;
  FILE *fid_in;
  float *map_in_f;
  double *map_in_d, *map_in_w;
  double dk, dk1, dk2, dk3, *pwsp3d;
  double complex *map_out;
  fftw_plan pr2c;
  int flag;
  double wsum;


  /* dl is now size of cell */
#ifdef _OMPTHREAD_
  printf("# SYNTAX: power3d Nbytes N1 N2 N3 dl flag input_file nthreads (output goes to standard output)\n\n");
  if(argc != 9) exit(1);
  nthreads=atoi(argv[8]);
  omp_set_num_threads(nthreads);
#else
  printf("# SYNTAX: power3d Nbytes N1 N2 N3 dl flag input_file (output goes to standard output)\n\n");
  if(argc != 8)  exit(1); 
#endif
  printf("#\"Nbytes\" is the number of bytes of each value (e.g., 4 for float and 8 for double).\n# The input file should have the following format: N3 is the fastest index, then N2, then N1 (e.g. the first set of numbers in the input binary file should be for N3 and so on...).\n# dl is the cell size (the same in all directions). The k and P(k) will have units after dl.\n# flag selects the window function: 0 - top hat, 1 - Bartlett, 2 - Welch\n");

#ifdef _OMPTHREAD_
  if(fftw_init_threads()==0) {
    printf("Error fftw threads!\n");
    exit(1);
  }
#endif
#ifdef _OMPTHREAD_
  fftw_plan_with_nthreads(nthreads);
#endif

  bsize=atoi(argv[1]);
  N1=(long int)atoi(argv[2]);
  N2=(long int)atoi(argv[3]);
  N3=(long int)atoi(argv[4]);
  dl=atof(argv[5]);
  flag=atoi(argv[6]);
  L1=N1*dl;
  L2=N2*dl;
  L3=N3*dl;
  dk1=2.0*PI/L1;  /* frequency size for power spectrum definition - k=2*PI*k' where k' is the FFT definition */
  dk2=2.0*PI/L2;
  dk3=2.0*PI/L3;
  /* take dk as the largest bin. In principle we could choose smaller dk or maybe a log dk instead... */
  if(dk1>dk2) dk=dk1; else dk=dk2;
  if(dk3>dk) dk=dk3; 

  printf("# Input parameters - Nbytes:%d, N1:%ld,  N2:%ld, N3:%ld, dl:%f, flag:%d, L1:%f, L2:%f, L3:%f, dk:%f, input file:%s\n", bsize, N1, N2, N3, dl, flag, L1, L2, L3, dk,argv[7]);
  
  nk=(long int)(round(sqrt(((N1/2)*(N1/2))+((N2/2)*(N2/2))+((N3/2)*(N3/2)))))+1;
  nps=(long int *)malloc(nk*sizeof(long int));
  pwsp3d=(double *)malloc(nk*sizeof(double));
  for(i=0;i<nk;i++) {
    nps[i]=0;
    pwsp3d[i]=0.;
  }

  fid_in=fopen(argv[7],"rb");	

  if(bsize==4) {
    if(!(map_in_f=(float *) malloc(N1*N2*N3*sizeof(float)))) {
      printf("Mem1 Problem...\n");
      exit(1);
    }
    fread(map_in_f,bsize,N1*N2*N3,fid_in);
  }
  if(bsize==8) {
    if(!(map_in_d=(double *) malloc(N1*N2*N3*sizeof(double)))) {
      printf("Mem1 Problem...\n");
      exit(1);
    }
    fread(map_in_d,bsize,N1*N2*N3,fid_in);
  }
  if(!(map_in_w=(double *) malloc(N1*N2*N3*sizeof(double)))) {
    printf("Memory Problem...\n");
    exit(1);
  }

  wsum=0.;
  for(i=0;i<N1;i++) {
    for(j=0;j<N2;j++) {
      for(p=0;p<N3;p++) {
	wsum+=wind(i,j,p,N1-1,N2-1,N3-1,flag)*wind(i,j,p,N1-1,N2-1,N3-1,flag);  /* normalize power spectrum (due to cut sky) */
      }
    }
  }
  printf("#sum: %E\n",wsum); fflush(0);

  if(bsize==4) { 
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_in_w,map_in_f,N1,N2,N3) private(i,j,p)
#endif
    for(i=0;i<N1;i++) {
      for(j=0;j<N2;j++) {
	for(p=0;p<N3;p++) {
	  map_in_w[i*N2*N3+j*N3+p]=(double)map_in_f[i*N2*N3+j*N3+p]*wind(i,j,p,N1-1,N2-1,N3-1,flag);  /* .....apply window to patch */
	}
      }
    }
    free(map_in_f);
  }
  if(bsize==8) { 
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(map_in_w,map_in_d,N1,N2,N3) private(i,j,p)
#endif
    for(i=0;i<N1;i++) {
      for(j=0;j<N2;j++) {
	for(p=0;p<N3;p++) {
	  map_in_w[i*N2*N3+j*N3+p]=map_in_d[i*N2*N3+j*N3+p]*wind(i,j,p,N1-1,N2-1,N3-1,flag);  /* .....apply window to patch */
	}
      }
    }
    free(map_in_d);
  }

  if(!(map_out = (double complex*) malloc(sizeof(double complex) * N1*N2*(N3/2+1)))) {
    printf("Mem3 Problem...\n");
    exit(1);
  }
  if(!(pr2c=fftw_plan_dft_r2c_3d(N1, N2, N3, map_in_w, map_out, FFTW_ESTIMATE))) {
    printf("FFT Problem...\n");
    exit(1);
  }
  fftw_execute(pr2c);

  printf("#Now calculating power spectum\n");fflush(0);
  /* Calculate power */
  for(i=0;i<N1;i++) {
    if(i>N1/2) {
      indi=-(N1-i);
    }else indi=i;
    for(j=0;j<N2;j++) {
      if(j>N2/2) {
	indj=-(N2-j);
      }else indj=j;
      for(p=0;p<=N3/2;p++) {	
	k=sqrt((indi*dk1)*(indi*dk1)+(indj*dk2)*(indj*dk2)+(p*dk3)*(p*dk3));
	indk=(long int)round(k/dk);
        nps[indk]++;
        pwsp3d[indk]+=(creal(map_out[i*N2*(N3/2+1)+j*(N3/2+1)+p])*creal(map_out[i*N2*(N3/2+1)+j*(N3/2+1)+p])+cimag(map_out[i*N2*(N3/2+1)+j*(N3/2+1)+p])*cimag(map_out[i*N2*(N3/2+1)+j*(N3/2+1)+p]));
      }
    }
  }

  for(i=0;i<nk;i++) {
    if(nps[i]==0) nps[i]=1;
      printf("%f   %E\n",i*dk,pwsp3d[i]/nps[i]/wsum*dl*dl*dl); 
  }

  fftw_destroy_plan(pr2c);
  free(map_out);
  free(map_in_w);
  fclose(fid_in);
  exit(0);

}


/* Window function */
double wind(int j,int k,int p,int Nx,int Ny,int Nz, int flag) {

  if(flag==1) return (1-fabs((j)-(Nx)/2.)/(Nx)*2.)*(1-fabs((k)-(Ny)/2.)/(Ny)*2.)*(1-fabs((p)-(Nz)/2.)/(Nz)*2.); // Bartlett 
  if(flag==2) return (1-SQR(((j)-(Nx)/2.)/(Nx)*2.))*(1-SQR(((k)-(Ny)/2.)/(Ny)*2.))*(1-SQR(((p)-(Nz)/2.)/(Nz)*2.)); // Welch 
  return 1.0; // top hat

}

