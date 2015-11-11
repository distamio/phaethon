#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "params.grid.h"
#include "params.h"
#include "constants.h"
#include "inout.h"
#include "opacity.h"
#include "phaethon.functions.h"
#include "grid.functions.h"


int cell_dustkind(ctypeptr, int);

/****** PHOTON DUST KIND************************************/

#ifdef DISC_IN_CLOUD

void find_photon_dust_kind(p,i)

     struct photon_packet_type p[];
     int i;
{

  if (DISC==1 && p[i].r<DISC_R) {p[i].dust_kind=0;}

  else    

{p[i].dust_kind=1;}



  return;
}

void find_cell_k_visual(c,i)

     struct cell_type c[];
     int i;
{

  if (DISC==1 && c[i].x<DISC_R) c[i].K_VISUAL=K_VISUAL0;

  else    

c[i].K_VISUAL=K_VISUAL1;

  return;
}

int cell_dustkind(c,i)

     struct cell_type c[];
     int i;
{

  if (DISC==1 && c[i].x<DISC_R) return(0);

  else    

 return (1);

  return;
}


double find_k_visual(r)

    double r;
{
 

  if (DISC==1 && r<DISC_R) return(K_VISUAL0);

  else    

  return(K_VISUAL1);
}


#endif
