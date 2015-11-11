/* OPACITY_SPH  opacity_OH5_0.1K_100K +  MRN/2 (Draine 03, not rescaled) */

/********************** pemit *****************************************/
#define NFREQ  1000  /* number of frequencies to reemit photons */

/* wavelength limits to reemit photon packets [pemit.c]*/

#define L1     0.1 /* microns */
#define L2     10000.


/*************************** opacity    ********************************/

#define MAX_TEMP  101.  /* max temperature to calculate mean opacity +1 */
#define TEMPDIM   1001    
  
#define OPADIM    1000 /* dim  >100 of the opacity.dat, if there is one    */
