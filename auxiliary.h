
/*
SimFast21
*/

#ifndef __AUXILIARY__
#define __AUXILIARY__

#include <complex.h>


/*-----------Type def for Halo information----------------------------------*/
typedef struct Halo_t_
{
  float Mass;
  float x,y,z;
} Halo_t;


/* For T_k */
typedef struct tf_parms_
{
  double omhh;
  double theta_cmb;
  double growth_k0;
  double alpha_gamma;
  double sound_horizon_fit;
  double beta_c;
  double growth_to_z0;
} tf_parms;


double W2(double x);

int Set_Cosmology(double omega_matter, double omega_baryon, double omega_lambda, double redshift, tf_parms *tf);
double getGrowth(double z);
double getGrowthb(double z);
double powerFunction(double k, tf_parms *tf);
double sig8(double omm, double omb, double lab);
double sigma(double R);
double deltaFilter(double sigma_aux, double growth);
int poisson_sampling(double dn, double growth, double sigma2, double deltarho, double z, double dx, double bias,double p);
double Bias(double z, double sigmaM);
double W_filter(double x);
double mass_function_ST(double z, double M);
long int check_borders(long int x,long int N);
double dzdt(double z);
double dgrowthdt(double redshift);
double dgrowthdz(double redshift);
double Hz(double z);


void smooth_box(float *box, float *box_smoothed, long int N, long int Ns);
void box_symmetriesd(double complex *box, long int N);
void box_symmetriesf(float complex *box, long int N);
void get_collapsed_mass_box(float* halo_box,Halo_t *halo, long int nhalos);
void get_collapsed_mass_boxb(float* halo_box,Halo_t *halo, long int nhalos);
void CIC(float x1, float y1, float z1, float map_in, float *map_out, long int N);



#endif


