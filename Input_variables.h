
/*
SimFast21
*/

#ifndef __INPUT_VARIABLES__
#define __INPUT_VARIABLES__


#define PI 3.1415926535897932385E0
#define RHO_0 2.775e11 /* units of Msun/h/(Mpc/h)^3 */


/*----Declaration of functions for reading params and setting vars-----------*/

void set_cosmology_fromCAMB(char * paramfilename); 
void get_Simfast21_params(char * basedir); 
void print_parms(void);


/*--------------------Simfast21 variables and parameters---------------------------*/

int global_nthreads;
long int global_seed;
long int global_N_halo; //Linear number of cells of the box for determination of collapsed halos 
long int global_N3_halo; //Total number of cells of the box for determination of collapsed halos 
long int global_N_smooth; // Linear number of cells of the smoothed boxes  
long int global_N3_smooth; // Total number of cells of the smoothed boxes  
float global_smooth_factor; //Just N_halo/N_smooth
double global_L; //Physical size of the simulation box
double global_L3;//Physical volume of the simulation box
double global_dx_halo;
double global_dx_smooth;
double global_dk;
int global_vi;
/* simulation range */
double global_Dzsim;
double global_Zmaxsim;
double global_Zminsim;
double global_Zminsfr;


/*-----------------------Cosmological parameters--------------------------- */

double global_sig8_new;
double global_n_index;
double global_hubble;
double global_omega_m;   
double global_omega_b;   
double global_lambda;
double global_rho_m;
double global_rho_b;

/*------------------------Halo collapse parameters-----------------------------*/

double global_delta_c;
double global_STa;
double global_STb;
double global_STc;
double global_halo_Rmax;
double global_halo_Rmin_dx;
double global_halo_Mmin;
int global_Nhbins;



/*------------------------Ionization parameters-----------------------------*/
double global_bubble_Rmax;
int global_bubble_Nbins;
double global_xHlim;   /* cutoff limit for bubble calculation */
double global_fesc; /* escape fraction */


/*----------Variables for reading matter power spectrum from file-------- */

int global_pk_flag; // Matter power spectrum: 0 - Eisenstein & Hu fitting formulae; 1 - Read form file (Output of CMBFast, CAMB)
char global_pk_filename[99];
char global_camb_file[99];



/*-------------------Flags for output files and algorithm----------------------------------*/ 
int global_use_sgrid;
int global_save_nl_halo_cat;
int global_save_original_deltanl;
int global_use_Lya_xrays;


/*--------- Parameters for X-ray heating and Lya coupling----------*/
double global_fstar;
double global_Enu0;
double global_alphas;
double global_L0; 		
double global_flux_Rmax;
double global_A_Lya;
double global_alpha_Lya;



/* Additional variable for Lyalpha coupling and X-ray heating */

#define G 6.67300e-11 
#define Msun 1.989e30 /* Kg */
#define mbar 1.67e-27 /* = mproton (Kg) */

#define mH (1.00794e-3/6.02214e23) /*  Kg */
#define mHe (4.002602e-3/6.02214e23) /* Kg */
#define YHe 0.24 /* Helium mass fraction - it might be a input of CAMB*/
#define y2s (365.25*24.*3600.) /* s */
#define Mpc2m 3.08568025e22 /* m */
#define H0 3.24078e-18 /* h/sec */
//#define nmax 16 /* put as a parameter ?*/


#endif

 

