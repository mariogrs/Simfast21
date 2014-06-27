
/*
SimFast21
Auxiliary functions (cosmology, etc)
*/

#ifdef _OMPTHREAD_
#include <omp.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "Input_variables.h"
#include "auxiliary.h"
#include <gsl/gsl_integration.h>

#define SQR(a) ((a)*(a))



/* ------------------------- TFmdm_set_cosm() ------------------------ */
int Set_Cosmology(double omega_matter, double omega_baryon, double omega_lambda, double redshift, tf_parms *tf) {
/* This routine takes cosmological parameters and a redshift and sets up
all the internal scalar quantities needed to compute the transfer function. */
/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. */
/* 	  omega_baryon -- Density of baryons, in units of critical. */
/*        omega_lambda -- Cosmological constant */
/* 	  global_hubble       -- hubble constant, in units of 100 km/s/Mpc */
/*        redshift     -- The redshift at which to evaluate */
/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
	sets many global variables for use in TFmdm_onek_mpc() */


/*------- Variables for Einsestein & Hu transfer function-----------------*/
  double alpha_gamma,     /* sqrt(alpha_nu) */
    alpha_nu,        /* The small*/
    beta_c,		   /* The correction to the log in the small-scale */
    f_baryon,	       /* Baryon fraction */
    f_cb,		       /* Baryon + CDM fraction */
    f_cdm,		   /* CDM fraction */
    growth_k0,	   /* D_1(z) -- the growth function as k->0 */
    growth_to_z0,	   /* D_1(z)/D_1(0) -- the growth relative to z=0 */
    k_equality,	   /* The comoving wave number of the horizon at equality*/
    obhh,		       /* omega_baryon * hubble^2 */
    omega_curv,	   /* = 1 - omega_matter - omega_lambda */
    omega_lambda_z,  /* omega_lambda at the given redshift */
    omega_matter_z,  /* omega_matter at the given redshift */
    omhh,		       /* omega_matter * hubble^2 */
    p_c,		       /* The correction to the exponent before drag epoch */
    p_cb,		       /* The correction to the exponent after drag epoch */
    sound_horizon_fit,  /* The sound horizon at the drag epoch */
    theta_cmb,	   /* The temperature of the CMB, in units of 2.7 K */
    y_drag,		   /* Ratio of z_equality to z_drag */
    z_drag,		   /* Redshift of the drag epoch */
    z_equality;	   /* Redshift of matter-radiation equality */
  

  /********************************************************************/
  
  double z_drag_b1, z_drag_b2, omega_denom;
  int qwarn;
  qwarn = 0;

  theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */
  
  /* Look for strange input */
  if (omega_baryon<0.0) {
    fprintf(stderr,
	    "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
    qwarn = 1;
  }
  
  if (global_hubble<=0.0) {
    fprintf(stderr,"TFmdm_set_cosm(): Negative hubble constant illegal.\n");
    exit(1);  /* Can't recover */
  } else if (global_hubble>2.0) {
    fprintf(stderr,"TFmdm_set_cosm(): hubble constant should be in units of 100 km/s/Mpc.\n");
    qwarn = 1;
  }
  if (redshift<=-1.0) {
    fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
    exit(1); //terminates the prog. the value 1 means that the program went well
  } else if (redshift>99.0) {
    fprintf(stderr,
	    "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
    qwarn = 1;
  }
  
  if (omega_baryon<=0) omega_baryon=1e-5;
  
  omega_curv = 1.0-omega_matter-omega_lambda;
  omhh = omega_matter*SQR(global_hubble);
  obhh = omega_baryon*SQR(global_hubble);
  f_baryon = omega_baryon/omega_matter;
  
  f_cdm = 1.0-f_baryon;      
  f_cb = 1;                   
  
  /* Compute the equality scale. */
  z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
  k_equality = 0.0746*omhh/SQR(theta_cmb);
  
  /* Compute the drag epoch and sound horizon */
  z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
  z_drag_b2 = 0.238*pow(omhh,0.223);
  z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*
    (1.0+z_drag_b1*pow(obhh,z_drag_b2));
  y_drag = z_equality/(1.0+z_drag);
  
  sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));
  
  /* Set up for the free-streaming & infall growth function */
  p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
  p_cb = 0;   
  
  omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+
						omega_matter*(1.0+redshift));     
  omega_lambda_z = omega_lambda/omega_denom;
  omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
  growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
    (pow(omega_matter_z,4.0/7.0)-omega_lambda_z+
     (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
  growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)
					      -omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
  growth_to_z0 = growth_k0/growth_to_z0;	
    
  /* Compute small-scale suppression */  
  alpha_nu = (1- f_baryon)* (5-2*p_c)/(5)* pow(1+y_drag,-p_c)
    *(1+f_baryon*(-0.553+0.126*f_baryon*f_baryon))*(1+(p_c)/2*(1+1/(3.-4.*p_c)
							       /(7.))/(1+y_drag));
  
  alpha_gamma = sqrt(alpha_nu); 
  beta_c = 1/(1-0.949*f_baryon);
  
  (*tf).omhh=omhh;
  (*tf).theta_cmb=theta_cmb;
  (*tf).growth_k0=growth_k0;
  (*tf).alpha_gamma=alpha_gamma;
  (*tf).sound_horizon_fit=sound_horizon_fit;
  (*tf).beta_c=beta_c;  
  (*tf).growth_to_z0=growth_to_z0;
    
  return qwarn;
}



/* ---------------------------- TFmdm_onek_mpc() ---------------------- */

double Transfer_Function(double kk, tf_parms *tf) {
  /* Given a wavenumber in Mpc^-1, return the transfer function for the
     cosmology held in the global variables. */
  /* Input: kk -- Wavenumber in Mpc^-1 */
  /* Output: the transfer function for density-weighted
     CDM + Baryon perturbations. */
  
  double tf_sup_L, tf_sup_C;
  double temp1, temp2;

  /* The following are set in TFmdm_onek_mpc() */
  double gamma_eff,	   /* Effective \Gamma */
    growth_cb,	   /* Growth factor for CDM+Baryon perturbations */
    //    max_fs_correction,  /* Correction near maximal free streaming */
    qq,		       /* Wavenumber rescaled by \Gamma */
    qq_eff,		   /* Wavenumber rescaled by effective Gamma */
    //    qq_nu,		   /* Wavenumber compared to maximal free streaming */
    //    tf_master,	   /* Master TF */
    tf_sup,		   /* Suppressed TF */
    y_freestream;    /* The epoch of free-streaming for a given scale */
  
  /* TFmdm_onek_mpc() give its answers as */
  //  double tf_cb;		  /* The transfer function for density-weighted CDM + Baryon perturbations. */

  double omhh,theta_cmb,growth_k0,alpha_gamma,sound_horizon_fit,beta_c;

  omhh=(*tf).omhh;
  theta_cmb=(*tf).theta_cmb;
  growth_k0=(*tf).growth_k0;
  alpha_gamma=(*tf).alpha_gamma;
  sound_horizon_fit=(*tf).sound_horizon_fit;
  beta_c=(*tf).beta_c;
  
  qq = kk/omhh*SQR(theta_cmb);
  
  /* Compute the scale-dependent growth functions */
  y_freestream = 0;
  temp1= growth_k0;
  temp2 = pow(growth_k0,0.7); 
  growth_cb = growth_k0;
  
  /* Compute the master function */
  gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
		   (1+SQR(SQR(kk*sound_horizon_fit*0.43))));
  qq_eff = qq*omhh/gamma_eff;
  
  tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
  tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
  tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));
    
  return tf_sup;
}


double Transfer_Functionh(double kk, tf_parms *tf){//Tk com k in units of hMpc-1 
      
  return Transfer_Function(kk*global_hubble,tf);
      
}



/* ---------------------------- Hz() ---------------------- */
/* units s^-1 */
double Hz(double z) {
  return global_hubble*H0*sqrt(global_omega_m*(1.+z)*(1.+z)*(1.+z)+global_lambda); 
}



/* dz/dt in 1/seconds */
double dzdt(double z){
  return (-Hz(z)*(1.0+z));
}



/* Returns the integral from z to z+100 */

double intGrowth(double zz){

  double dz=0.001;
  double intg=0.;
  double Hubz,z,Hz0;

  Hz0=Hz(0.);
  for(z=zz+dz/2.;z<zz+100.;z+=dz) {
    Hubz=Hz(z)/Hz0;
    intg+=(1.+z)/Hubz/Hubz/Hubz;
  }
  intg=intg*dz;

  return intg;
}



double getGrowth(double z){

  //  printf("%E %E int: %E  %E\n",Hz(z),Hz(0),intGrowth(z),intGrowth(0));
  return Hz(z)/Hz(0)*intGrowth(z)/intGrowth(0.);
}

double getGrowthb(double z){

  int   eSetCosm;          // = 0 if their is an error = 1 if  otherwise
  tf_parms tf;

  /*Inicializacao dos parametros cosmologicos*/ 
  eSetCosm= Set_Cosmology(global_omega_m, global_omega_b, global_lambda, z, &tf);
    
  return tf.growth_to_z0;
}


double dgrowthdt(double redshift){
  double dg;
  double dz= 1e-3;

  dg = (getGrowth(redshift+dz/2.)-getGrowth(redshift-dz/2.))/dz*dzdt(redshift);
 
  return dg;
}

double dgrowthdz(double redshift){
  double dg;
  double dz= 1e-3;

  dg = (getGrowth(redshift+dz/2.)-getGrowth(redshift-dz/2.))/dz;
 
  return dg;
}




/*----------------------Função Power----------------------*/
/* P(k) for z=0 */
/* Assumes Set_Cosmology() was called */
double powerFunction(double k, tf_parms *tf){
  
  double tfk,        // Transfer function result
    power;         // Power spectra result
  double deltaH=4.2E-5; /* Cobe Norm */
  
  //k in units of hMpc-1 
  tfk = Transfer_Functionh(k,tf);
  power = (*tf).growth_to_z0*(*tf).growth_to_z0*(pow((1/k),3)*2*pow(PI*deltaH,2)*pow((k*global_hubble),global_n_index+3)*pow(tfk,2));
  //power in units of de (Mpc/h)^3
  
  return power;
}



/*--------------------------------------------------------------------------------*/
double W2(double x){ 
  
  double j1onx;
  
  if (x<0.03) 
    j1onx = (1.0/3.0-x*x*1.0/30.0);
  else 
    j1onx = -(cos(x)-(sin(x)/x))/x/x;
  return 9.0*j1onx*j1onx;  //9*j1(kr)/KR*j1(kr)/kR=W2
}

double W_filter(double x){ 
  
  double j1onx;
  
  if (x<0.03) 
    j1onx = (1.0/3.0-x*x*1.0/30.0);
  else 
    j1onx = -(cos(x)-(sin(x)/x))/x/x;
  return 3.0*j1onx;  
}


double dsigmaR(double x, void *parms)
{ 
  double R=*(double *)(*(void **)parms);
  tf_parms *tf=(tf_parms *)(*((void **)parms+1));

  return (powerFunction(x/R,tf)/(2*PI*PI)*(pow(x/R,3))*W2(x)/x+
	  powerFunction(1/x/R,tf)/(2*PI*PI)*(pow(1/x/R,3))*W2(1/x)/x);
}


/* R in Mpc/h */
double sigmaR(double R, double omm, double omb, double lamb) {
  double in_result,in_error;
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  int eSetCosm;
  tf_parms tf;
  void *parms[2];

  parms[0]=&R;
  parms[1]=&tf;

  eSetCosm= Set_Cosmology(omm, omb, lamb, 0.0, &tf); /* sets cosmology for transfer function calculation at z=0 */

  F.function = &dsigmaR;
  F.params = parms;

  gsl_integration_qag (&F, 0., 1.,0., 3e-7, 1000,
                       GSL_INTEG_GAUSS15,
                       w, &in_result, &in_error);
  gsl_integration_workspace_free (w);
  return (sqrt(in_result));
}


/* sigma8 with Cobe normalization */
double sig8(double omm, double omb, double lamb) {

  return sigmaR(8.0,omm,omb,lamb);  /* 8 Mpc/h */

}
  

double sigma(double R) {

  double sig8_old;

  sig8_old=sig8(global_omega_m, global_omega_b, global_lambda);

  return global_sig8_new/sig8_old*sigmaR(R,global_omega_m, global_omega_b, global_lambda);

}


double deltaFilter(double sigma_aux, double growth){

  double deltacMZ;
  double dz;
  
   
  dz=global_delta_c/growth;
  
  deltacMZ=sqrt(global_STa)*dz*(1.0+global_STb*pow((sigma_aux*sigma_aux/(global_STa*dz*dz)),global_STc));
        
  return deltacMZ;    
}



double Bias(double z, double sigmaM){

  double b,deltac,growth,nuhat;

  deltac =global_delta_c;
  growth=getGrowth(z);
  sigmaM=sigmaM*growth;
  nuhat=deltac/sigmaM;
  b=(1.0+((nuhat*nuhat-1.0)/deltac))*pow((1.0+(1/(2.0*pow(nuhat,4)))),0.06-0.02*global_n_index);

  return b;
}


double poisson_mean(double dn,double growth, double sigma2, double deltarho, double z,double dx,double bias){

  double poisson,growth2,B;  
 
  growth2=growth*growth;
  //A=exp(-bias*bias*growth2*sigma2/2.);
  B=-bias*bias*growth2*sigma2/2. + bias*deltarho;
  
  poisson=exp(B)*dn;
  //poisson=dn;
  //printf("poisson_mean=%lf\n",poisson);fflush(0);
  return poisson;
}


int poisson_sampling(double dn, double growth, double sigma2, double deltarho, double z, double dx, double bias,double p){
  
  int check;
  double mean;

  mean=poisson_mean(dn,growth,sigma2,deltarho,z,dx,bias);
  if(mean>=p)check=1;
  else check=0;
  return check;
  
}


/* output in Mpc/h */
double MasstoRadius(double M){
  return pow(3*M/(4*PI*global_omega_m*RHO_0/global_hubble), 1.0/3.0);
}


/* Mass in Msun */
double mass_function_ST(double z, double M){

  double deltac,sigmaM,sigmaMdm, dlnsigmadm, nuhat, growth,M2,R1,R2;
  double A,a,p, result,E;

  M2=M+(M/1000.); 
  a=0.73;
  A=0.353;
  p=0.175;
  E=2.71828182846;
  R1=MasstoRadius(M);
  R2=MasstoRadius(M2);
    
  growth = getGrowth(z);

  deltac=global_delta_c;
 
  sigmaM = sigma(R1)*growth ; 
  sigmaMdm = sigma(R2)*growth; 
  
  dlnsigmadm = (sigmaMdm-sigmaM)/(M/1000.);
 
  nuhat = sqrt(a) * deltac / sigmaM;
  
  result=(-global_omega_m*RHO_0/global_hubble)/M * dlnsigmadm/sigmaM * sqrt(2.0/PI)*A * (1+ pow(nuhat, -2*p)) * nuhat * pow(E, -nuhat*nuhat/2.0);
 
  return result;

}



long int check_borders(long int x, long int N){
  
  x=x%N;
  if(x<0)x+=N;

  return x;
  
}



/* ----------------Functions---------------------------------------*/
void box_symmetriesd(double complex *box, long int N) {

  long int i,j,indNi,indNj;

  for(i=0;i<N;i++) {  // do box[i,j,0]
    for(j=0;j<N;j++) {
      if(i==0) {
	indNi=0;
      }else indNi=N-i;
      if(j==0) {
	indNj=0;
      }else indNj=N-j;
      box[indNi*N*(N/2+1)+indNj*(N/2+1)+0]=conj(box[i*N*(N/2+1)+j*(N/2+1)+0]);
    }
  }
  box[0]=0.;
  if(N%2==0) {    // N is even - do also box[i,j,N/2]
    for(i=0;i<N;i++) {
      for(j=0;j<N;j++) {
	if(i==0) {
	  indNi=0;
	}else indNi=N-i;
	if(j==0) {
	  indNj=0;
	}else indNj=N-j;
	box[indNi*N*(N/2+1)+indNj*(N/2+1)+N/2]=conj(box[i*N*(N/2+1)+j*(N/2+1)+N/2]);
      }
    }
    i=N/2;
    box[i*N*(N/2+1)+0*(N/2+1)+0]=creal(box[i*N*(N/2+1)+0*(N/2+1)+0])+I*0.;
    box[0*N*(N/2+1)+i*(N/2+1)+0]=creal(box[0*N*(N/2+1)+i*(N/2+1)+0])+I*0.;
    box[0*N*(N/2+1)+0*(N/2+1)+i]=creal(box[0*N*(N/2+1)+0*(N/2+1)+i])+I*0.;
    box[i*N*(N/2+1)+i*(N/2+1)+0]=creal(box[i*N*(N/2+1)+i*(N/2+1)+0])+I*0.;
    box[i*N*(N/2+1)+0*(N/2+1)+i]=creal(box[i*N*(N/2+1)+0*(N/2+1)+i])+I*0.;
    box[0*N*(N/2+1)+i*(N/2+1)+i]=creal(box[0*N*(N/2+1)+i*(N/2+1)+i])+I*0.;
    box[i*N*(N/2+1)+i*(N/2+1)+i]=creal(box[i*N*(N/2+1)+i*(N/2+1)+i])+I*0.;
  }

}




void box_symmetriesf(float complex *box, long int N) {

  long int i,j,indNi,indNj;

  for(i=0;i<N;i++) {  // do box[i,j,0]
    for(j=0;j<N;j++) {
      if(i==0) {
	indNi=0;
      }else indNi=N-i;
      if(j==0) {
	indNj=0;
      }else indNj=N-j;
      box[indNi*N*(N/2+1)+indNj*(N/2+1)+0]=conjf(box[i*N*(N/2+1)+j*(N/2+1)+0]);
    }
  }
  box[0]=0.;
  if(N%2==0) {    // N is even - do also box[i,j,N/2]
    for(i=0;i<N;i++) {
      for(j=0;j<N;j++) {
	if(i==0) {
	  indNi=0;
	}else indNi=N-i;
	if(j==0) {
	  indNj=0;
	}else indNj=N-j;
	box[indNi*N*(N/2+1)+indNj*(N/2+1)+N/2]=conjf(box[i*N*(N/2+1)+j*(N/2+1)+N/2]);
      }
    }
    i=N/2;
    box[i*N*(N/2+1)+0*(N/2+1)+0]=crealf(box[i*N*(N/2+1)+0*(N/2+1)+0])+I*0.;
    box[0*N*(N/2+1)+i*(N/2+1)+0]=crealf(box[0*N*(N/2+1)+i*(N/2+1)+0])+I*0.;
    box[0*N*(N/2+1)+0*(N/2+1)+i]=crealf(box[0*N*(N/2+1)+0*(N/2+1)+i])+I*0.;
    box[i*N*(N/2+1)+i*(N/2+1)+0]=crealf(box[i*N*(N/2+1)+i*(N/2+1)+0])+I*0.;
    box[i*N*(N/2+1)+0*(N/2+1)+i]=crealf(box[i*N*(N/2+1)+0*(N/2+1)+i])+I*0.;
    box[0*N*(N/2+1)+i*(N/2+1)+i]=crealf(box[0*N*(N/2+1)+i*(N/2+1)+i])+I*0.;
    box[i*N*(N/2+1)+i*(N/2+1)+i]=crealf(box[i*N*(N/2+1)+i*(N/2+1)+i])+I*0.;
  }

}


float *smooth_boxb(float *box, float *box_smoothed, long int N, long int Ns) {

  double av,sm,sm3;
  long int i,j,p,ii,jj,pp,indi,indj,indp,ism;
  long int inds;

  if(Ns>N) {
    printf("Problem - smooth_box: output N larger than input N\n");fflush(0);
    exit(1);
  }

  sm=(1.0*N)/Ns;  /* Try an approximate thing in case N is not a multiple of Ns... */
  ism=(long int)sm;
  sm3=1.0*ism*ism*ism;
  //  printf("smoothing box: %ld %ld %ld %lf %lf\n",N,Ns,ism,sm,sm3);fflush(0);
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(N,Ns,box,box_smoothed,sm,ism,sm3) private(i,j,p,inds,indi,indj,indp,av,ii,jj,pp)
#endif 
  for(i=0;i<Ns;i++) {
    for(j=0;j<Ns;j++) {
      for(p=0;p<Ns;p++) {
	inds=i*Ns*Ns+j*Ns+p;
	indi=(long int)(i*sm);
	indj=(long int)(j*sm);
	indp=(long int)(p*sm);
	av=0.0;
	for(ii=0;ii<ism;ii++) {
	  for(jj=0;jj<ism;jj++) {
	    for(pp=0;pp<ism;pp++) {
	      av+=box[(indi+ii)*N*N+(indj+jj)*N+(indp+pp)];
	    }
	  }
	}
	av=av/sm3;
	box_smoothed[inds]=(float)av;
      }
    }
  }

  return box_smoothed;

}



void get_collapsed_mass_box(float* halo_box,Halo_t *halo, long int nhalos){
  
  long int ii,ij,ik, ii_c,ij_c, ik_c, a, b, c;
  long int il,ncells_1D, ncells_3D;

  /* Need to set halo_box to zero!!! */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_smooth,halo_box) private(ii)
#endif 
  for (ii=0;ii<global_N3_smooth;ii++) halo_box[ii]=0.;

  for(il=0;il<nhalos;il++) {
   
    ii_c=(long int)((halo[il].x)/global_smooth_factor); 
    ij_c=(long int)((halo[il].y)/global_smooth_factor);
    ik_c=(long int)((halo[il].z)/global_smooth_factor);
  
    ncells_1D=(long int)(halo[il].Radius/global_dx_smooth);
    if(ncells_1D==0)ncells_1D=1;
        
    ncells_3D=0;
    for(ii=-(ncells_1D+1);ii<=ncells_1D+1;ii++){
      for(ij=-(ncells_1D+1);ij<=ncells_1D+1;ij++){
	for(ik=-(ncells_1D+1);ik<=ncells_1D+1;ik++){	
	  if((ii*ii+ij*ij+ik*ik)*global_dx_smooth*global_dx_smooth <= halo[il].Radius*halo[il].Radius)ncells_3D++;
	}
      }
    }
    
    for(ii=-(ncells_1D+1);ii<=ncells_1D+1;ii++){
      a=ii_c+ii;
      a=check_borders(a,global_N_smooth);
      for(ij=-(ncells_1D+1);ij<=ncells_1D+1;ij++){
	b=ij_c+ij;
	b=check_borders(b,global_N_smooth);
	for(ik=-(ncells_1D+1);ik<=ncells_1D+1;ik++){
	  c=ik_c+ik;
	  c=check_borders(c,global_N_smooth);
	  if((ii*ii+ij*ij+ik*ik)*global_dx_smooth*global_dx_smooth <= halo[il].Radius*halo[il].Radius){
	    halo_box[a*global_N_smooth*global_N_smooth+b*global_N_smooth+c]+=halo[il].Mass/ncells_3D;  /* After Halos are adjusted there could be some overlap?? */	    
	  }
	}
      }
    }

  }

}



void get_collapsed_mass_boxb(float* halo_box,Halo_t *halo, long int nhalos){
  
  long int i,j,p;
  long int il;

  /* Need to set halo_box to zero!!! */
#ifdef _OMPTHREAD_
#pragma omp parallel for shared(global_N3_smooth,halo_box) private(il)
#endif 
  for (il=0;il<global_N3_smooth;il++) halo_box[il]=0.;

  for(il=0;il<nhalos;il++) {    
    i=(long int)((halo[il].x)/global_smooth_factor); 
    j=(long int)((halo[il].y)/global_smooth_factor);
    p=(long int)((halo[il].z)/global_smooth_factor);
    halo_box[i*global_N_smooth*global_N_smooth+j*global_N_smooth+p]+=halo[il].Mass;  /* After Halos are adjusted there could be some overlap?? */	
  }

}    




