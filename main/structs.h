#include "params.h"
#include "phaethon.params.h"
#include "params.opa.h"
#include "params.grid.h"
#include "constants.h"

#define FPOINTS  NFREQ
#define MAX_ID 4
#define ISRF_MAX_DIM 80000

#ifdef BE_SPHERE
struct ftau_type
{
 double x[FTAU_POINTS];
 double tau[FTAU_POINTS];
};
#endif


struct star_type                      /* star (+disc) parameters or isrf parameters */
                                      /* could be star or isrf depending on star[].type */
            {
  
             double x;  
             double y;
             double z; 
	     unsigned short  sg;                 /* yes or no star grid */
	     unsigned short type;         /*0: star in the centre of the core
						    			1: isrf (Black or...) 
						   			    2: diluted star isrf 
						    			3: star outside the cloud (on Z axis)
						    			4: isrf 1D  
						 */
	      
	     double r_dust;           /* dust destruction radius */
	     double r_star_grid;      /* radious of star grid    */
	     double disc_kx;  
	     double disc_ky;
	     double disc_kz; 
	     unsigned long firstcellid;          /* id of the first cell of the star grid  */
	     unsigned long lastcellid;           /* id of the last  cell of the star grid  */

	     double R;                /* distance from (0,0,0) */
	     double mass;
             double temp;  
             double radius;        
   
             double   npackets;       /* number of photons emitted from this star */
	     unsigned long   cell;             /* cell where the star is in */
	     unsigned  short disc;                /* 1: star with disc 0: no disc */
	     double disc_z;
	     double disc_R;

	     double disc_mass;
	     double lum;               /*star luminosity */
	      unsigned short  rcells;               /* number of radial cells of the sphere/disc */
	      unsigned short lcells;               /* number of l cells of the disc */
	      unsigned short mcells;               /* number of m cells (just for postrt) */
            
	     double disc_theta_opening;/* opening angle from the disc midplane */
	     double disc_alpha;        /* z/R for disc */
	     double dilution;          /* dilution factor for the star radiation */

	      double output;               /* when to output info */
	      double tint;             /* interactions table: max interactions 
					             5000 for star 50 for isrf */
	      int dint;                /* interval 100 star, 1 isrf */

	      double vx;               /* velocity */
	      double vy;
	      double vz;

              double sphtime;           /* time of the sph snapshot */
	      long sphstep;           /* sph step                 */
	      long sphid;             /* first sink=1 second sink =2 ... */
	      long sphindex;          /* index of the sph particle that has become sink */
              double mdot;            /* accretion rate */
	                              /* table with information about the SED */
	      double lambda[ISRF_MAX_DIM+1];     /* microns */
	      double intensity[ISRF_MAX_DIM+1];     
	      int    sed_dim;                      /* dimenstion of the SED table */
	      double int_total;                /* total intensity of the radiation field */
              unsigned short    sedfile; /* 1: you read from file the SED */
      };

typedef struct star_type *startypeptr;

struct photon_packet_type
           {
             double f;            /* photon frequency */
             double x;            /* photon position  */
             double y;
             double z; 
	     double dist_to_star; /* distance from the star */
             double r;            /* distance from the origin of the system */       
             double kx;           /* photon direction k: unit vector */
             double ky;   
             double kz;           
             double k_abs;        /* absorption opacity (depends just on the f) */
             double s_scat;
	     double albedo;
             double tau;          /* optical depth */
             double lum;          /* photon luminosity */
             double starlum;      /* total luminosity of the source star/radiation field */
             long   cell;         /* index of the cell that the photon is in */
	                          /*-1 if the photon is outside the grid */
             unsigned  long abs;  /* how many times the photon is absorbed */
             unsigned  long scat; /* how many times the photon is scattered */
             unsigned  short id;  /* photon id */
	                      	    /* 
			   	     id=0 all photons
			   	
				     id=1 direct source photons
				     id=2 flagged photons 
				     id=3 scattered
                                     id=4 reprocessed radiation
			          */
	     double theta;
	     double phi;             
             unsigned short flag;
	     unsigned short origin;/* 0 background 1 source */
	     unsigned int fpoint;  /* index of the frequency bin the photon is in */
             long    lcell;        /* last cell before the photon escapes */
             double lr;            /* last r before the photon escapes */
             double lx; 
             double ly; 
             double lz; 
             unsigned short    n;  /* freq bin maybe just have only this and not f for the photons ! */
             unsigned short    l;  /* theta bin */
	     unsigned short    m;  /* phi bin */
             double    multiplicity;
	     double empty_cell_size;
#ifdef SPHERICAL
             double b;             /* impact factor */  
#endif
             double weight;        /* photon weight */
             double ds;
	     double intensity;	     


	     unsigned short dust_kind;/* the kind of dust that the photon interacts 0, 1 */
#ifdef DUST_MIX
	     double k_abs2;      /* second absorption opacity (depends just on the f) */
             double s_scat2;     /* for dust mixtures  */
             double albedo2;
#endif   
	     double temp; /* temp of the cell reemiting the photon */
           };   


typedef struct photon_packet_type *phtypeptr;

struct opacity_type
       {
        double lambda[OPADIM];
        double abs[OPADIM];
        double scat[OPADIM];
        double albedo[OPADIM];
       };


struct opacity_table
       {
        double lambda[FPOINTS];
        double abs[FPOINTS];
        double scat[FPOINTS];
        double albedo[FPOINTS];
       };

/* sph particles */

struct particle_type
           {
             double x;  
             double y;
             double z; 
             double vx; 
             double vy;
             double vz;
             double rho;
             double mass;
             double temp;
             int iwas;
             double h; 
             long index;
	     double r;
             double lum;
#ifndef LOWMEM
	     double penergy;
	     double column2;
	     double column;
	     int    rbin; 	     
             double senergy;
	     double mu;
	     double mopa;
	     double mopap;
	     double atemp;
	     double nctime;
	     double rcool;
	     double rheat;
             int    step;
             double ctime;
             double htime;
             double sigma;
	     double qtoom;              
             int wbin; 
	     int zbin; 
	     long id;
	     double xx;
	     double yy;
	     double z1;
	     double z2;
#endif
             };
typedef struct particle_type *ptypeptr;


struct cell_type
     {
       double x;              /* x coordinate of the center of the cell */
       double y;              /* these could also be r theta phi ... */
       double z;
       double dens;           /* cell density */
 
       double nsph;           /* number of sph particles that are inside the cell */
       unsigned long  nscat;  /* total number of photon packets scattered */
       
       double temp;           /* temperature of the cell */
       double abslum;         /* total luminosity of the absorbed packets */
       unsigned long    nabs; /* total number of photon packets absorbed */
       double mass;           /* cell mass */
       long   next_n;         /* absorbing next_n photons will raise the temp to the next value */

       double tau;
       double size;           /* cell size */
       unsigned short newcell;/* 0: no sub cells in the cell 1: yes */
       long   subcell;        /* the index of the first subcell in the cell */
       double advsize;        /* ADVCACE cell size */
#ifdef POST_RT
       double prt_tau[NLBINS][NMBINS][8]; /* tau (visual) from cell till the end of system at different viewing angles (8 per cell)*/
#endif
       unsigned short    flag;   /* 0 before the first absorption, 1 after that */
       unsigned short  K_VISUAL;  /* the visual extintion opacity kind of dust in cell */

};

typedef struct cell_type *ctypeptr;


/* structure with basic parameters of a grid */

struct grid_type
       
       {
        double xmin;    /* min max values of the grid */
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        double dx;      /* width of the grid cells */
        double dy;
        double dz;
        int    xcells;  /* number of cells in the x direction */
        int    ycells;
        int    zcells;
        int    ncells;
        };


/* this structure contains the mean opacity table */

struct mean_op_type
       {
        double temp[TEMPDIM];
        double opacity[TEMPDIM];
       };


/* reemit table with temp, freq and probabilities */

struct reemit_type

{
 double temp[TEMPDIM];
 double freq[NFREQ];
 double prob[TEMPDIM][NFREQ];
 int last[TEMPDIM];
 double rfreq[TEMPDIM]; 
 int fpoint[TEMPDIM]; 
};

struct xtable_type
       {
       double freq[FPOINTS]; /* frequencies to emit the photon packets */
       long num[FPOINTS];     /* number of particles at each frequency */ 
       int multiplicity[FPOINTS]; 

       };



struct btype
       {
	 double lambda[ISRF_MAX_DIM+1];     /* microns */
	 double intensity[ISRF_MAX_DIM+1];     
	 int    dim; /* dimenstion of the iSRF table */
	 double int_total; /* total intensity of the radiation field */
       };

/* frequency bins */

// --------------
struct fbin_type
       {
	 double lambda1[NFBINS+1];
	 double lambda2[NFBINS+1];
	 double lambda[NFBINS+1];
	 /* wavelenth in microns */
	 double lum[NFBINS+1][MAX_ID+1][NLBINS][NMBINS];    /* luminosity of the photons in the bin */
	 double npackets[NFBINS+1][MAX_ID+1][NLBINS][NMBINS];  /* number of photon packets in the bin */
	 double opacity[NFBINS+1]; /* total opacity at the bin's frequency */ 
	 double abs[NFBINS+1];
	 double scat[NFBINS+1];
	
};                           

/* impact parambeter  bins */


struct bbin_type
{
         double b[NBBINS+1];         /* impact factor (cgs) */
	 double lum[NBBINS+1][NFBINS+1][MAX_ID+1];       /* luminosity of the photons in the bin */
	 double npackets[NBBINS+1][NFBINS+1][MAX_ID+1];     /* number of photon packets in the bin */
         double intensity[NBBINS+1][NFBINS+1][MAX_ID+1];
         double sigma[NBBINS+1]; /* surface density */
       };


/* direction direction bins */


struct dbin_type
{
  double theta[NLBINS][NMBINS];
  double phi[NLBINS][NMBINS]; 
  double domega[NLBINS][NMBINS]; 
};



struct isogrid_type
       {
	 double x[NXPIXELS][NYPIXELS];         
	 double y[NXPIXELS][NYPIXELS];      
	 double npackets[NXPIXELS][NYPIXELS][FS][MAX_ID+1][NLBINS][NMBINS];
         double intensity[NXPIXELS][NYPIXELS][FS][MAX_ID+1][NLBINS][NMBINS];
	 double tau[NXPIXELS][NYPIXELS][NLBINS][NMBINS]; 
	 double dx; /* pixel size */
	 double dy; /* pixel size */
       };

//-------------

#ifdef BE_SPHERE

#define BEPOINTS    20000   /* 42000 */  

struct besphere_type

{
  double x[BEPOINTS+2];    /* ksi        */
  double y[BEPOINTS+2];    /* psi        */
  double r[BEPOINTS+2];    /* radius cgs */
  double rho[BEPOINTS+2];  /* cgs g/cm^3 */
  double mass[BEPOINTS+2]; /* cgs        */
  double p_ksi[BEPOINTS+2];
  double mass_ksi[BEPOINTS+2];
  double R_ksi[BEPOINTS+2];
  double temp; 
  double rout;  /* AU         */
  double cdens; /* cgs g/cm^3 */
  double mass_tot; /* Msun       */ 
  double pext;  /* P4         */
  double mu;
  double ksi_out;
};

#endif


