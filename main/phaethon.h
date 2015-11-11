#include "params.h"
#include "phaethon.params.h"
#define MAX_ID  4  

double lambda_obs[FS]; 

/* wavelengths to construct isophotal maps */

double dlambda_obs[FS];
 
struct btype bisrf;
struct photon_packet_type phtype;
struct star_type sttype;

struct cell_type ctype;
struct opacity_type opacity;
struct mean_op_type mean_op;
struct reemit_type reemit;

struct opacity_table opa;

#ifdef DUST_MIX
struct opacity_type opacity2;
struct mean_op_type mean_op2;
struct reemit_type reemit2;
struct opacity_table opa2;
#endif

struct dbin_type dbin;


/* long   idum=-1;*/         /* initial value for random number generator                */
double mfp_free=0;      /* distance the photons travel in a step in  empty space    */
long    ncells=0;        /* number of cells in the grid -1 to run from 0 to NCELLS-1 */
double interactions=0.;  /* total number of interactions                             */
double time0=0.;         /* starting cpu clock time */
double time1=0.;         /* ending cpu clock time   */
double time2=0.;
double intercount[51];  /* percentage of photons that interact 0,1,2...,20 times    */
double source_lum=0.;   /* luminosity of the source */

struct fbin_type fbin;  /* frequency bins */

#ifdef SPHERICAL
struct bbin_type bbin;  /* impact factor bins (only for spherical case) */
#endif


struct isogrid_type isogrid; /* isophotal map grid */


double    total_fbin=0.;      /* number of photons placed in frequency bins     */
double    total_bbin[10];     /* number of photons placed in b bins             */
int       fbin_obs[10];  /* contains the fbin index of the considered wavelengths */
int       fbin_obs_max[10];  
int       fbin_obs_min[10];  

double    theta_observer;
double    phi_observer;

double    background[10];

double    domega[NLBINS][NMBINS];

/* double  dx=(XMAX-XMIN)/(1.*NXPIXELS);
   double  dy=(YMAX-YMIN)/(1.*NYPIXELS); */

unsigned short s;                          /* luminosity source counter */
int k_fpoints=0;               /* number of photon in xtable: needs to be zeroed for each lum source */
double output_counter=0;
double photon_i=1;                /* photon counter             */  
double rtime0=0;                /* time of the restart file */
int s_start=0;                  /* id of the first luminosity source 
				   (=0 for a normal run, !=0 for a restart) */
double output_counter_start=0;  /* (=0 for a normal run,!=0 for a restart) */
long k_fpoints_start=0;         /* (=0 for a normal run,!=0 for a restart) */
double photon_i_start=0;     
long xnum=0;
int sed_flag=0;                 /* */
int iso_flag=0;                 /* */
int pk_fpoint=0;
long pk_num=0;

double dragon_time;             /* time of the sph output file */
long dragon_step;
short star_flag=0;      /* 0: no stars in simulation */
