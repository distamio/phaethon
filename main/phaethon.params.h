#include "constants.h"
#include "params.grid.h"      

/***** GENERAL PAPAMETERS   *****************/

#define MAX_CELLS           600000
#define MAX_PHOTON_PACKETS  100
#define TRUE  1
#define FALSE 0

/********** POST SPH treatment params  ***********************************/

#define DO_TAU_ONLY   0        /* needs to be 0 when doing others */
    
#define POST_SED      1
#define DO_STAR_CELLS 1
#define DO_SPH_CELLS  1
#define DO_STAR       1 
#define POST_ISO      0

/***** SPH-RT PARAMETERS **************************************************/

/* #define SPH_TREE */

#define STAR_GRID     1
#define STAR_DISC     1

#define SPH_PARTICLES_PER_CELL 50


/****  SPH DISC PARAMS    ************************************************/
 
#define F                4. /* scale length factor */
#define STAR_GRID_FACTOR 1. /* determines how big is the star grid */
#define STAR_RCELLS  18   
#define STAR_LCELLS  31     /* >= 5 should be odd number */ 
#define STAR_MCELLS  30     /* just for post-RT analysis */
  


#define ORIGIN_SIZE (5500*AU)
/* #define SUB_REGION_SIZE   (12*AU) */ /* defined it acts as boundary of the system */
/* #define SUB_REGION_SIZE_X (210*AU)*/ 
/* #define SUB_REGION_SIZE_Y (210*AU)*/
/* #define SUB_REGION_SIZE_Z (210*AU)*/

/* #define SHADOW_TAU_MAX 20    */ 

