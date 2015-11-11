#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "params.grid.h"
#include "params.opa.h"
#include "constants.h"
#include "params.h"
#include "phaethon.params.h"

#include "functions.h"
#include "inout.h"
#include "opacity.h"
#include "phaethon.functions.h"
#include "photons.h"
#include "spec.functions.h"


/**************** make step info file    ********************************/

void make_step_info(time0_ptr,time1_ptr)

     double *time0_ptr;
     double *time1_ptr;


{  
  FILE *fp;
  
  *time0_ptr=time(NULL);
  *time1_ptr=time(NULL);  

  if (RESTART==1)
    {
      if (( fp = fopen("step.info", "a+")) == NULL)
	{
	  fprintf(stderr, "FAILED TO OPEN step.info FILE. \n");
	  exit(0);
	}
      
      fclose(fp);
    }
  else
  {
      if (( fp = fopen("step.info", "w")) == NULL)
	{
	  fprintf(stderr, "FAILED TO OPEN step.info FILE. \n");
	  exit(0);
	}
      
      fclose(fp);
    }

      return;
}

/**************** update interaction table ********************************/

void  update_interaction_table(p,i,intercount_ptr,TINT,DINT)
  
     struct photon_packet_type p[];
     long i;
     double *intercount_ptr;
     int TINT;
     int DINT;
     
{
  int t;
  int k=0;
  
  for(t=0;t<TINT;t=t+DINT) 
    {
      if ( ((p[i].abs+p[i].scat)>=t) && ((p[i].abs+p[i].scat)<(t+DINT))) 
	{
	  *(intercount_ptr+k)=*(intercount_ptr+k)+1;	

	  break;   
	} 
      k=k+1;	     
    }
  
  if (p[i].abs+p[i].scat>=TINT) (*(intercount_ptr+50))=(*(intercount_ptr+50))+1;
  
  return;
}

/**************** record interaction table ********************************/

void record_interactions(intercount_ptr,npackets,interactions,s,DINT,npackets_current)

     double *intercount_ptr;
     double npackets;
     double interactions;
     int s;
     int DINT;
     long npackets_current;
{
  int t;char filename[30];
  FILE *fp;
  
  sprintf(filename,"intercount.%d.dat",s);

  if (( fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN intercount.dat FILE. \n");
      exit(0);
    }
  fprintf(fp,"@photons DONE: %ld\n",npackets_current);
  fprintf(fp,"@Total interactions: %.0f\n",interactions);
  fprintf(fp,"@Interactions   percent of photons\n");
  
  for(t=0;t<=50;t++)  fprintf(fp,"  %5d    % 9.5f\n",t*DINT, 100.*(*(intercount_ptr+t))/npackets); 
  
  fclose(fp);
  
  return;
}



/********* update step info file **********************************************/

void update_step_info(photon_i,npackets,time0, time1_ptr,s,rtime0)
     double photon_i;
     double npackets;
     double time0;
     double *time1_ptr;
     int s;
     double rtime0; /* restart end time */

{
  FILE *fp;
  
  if (( fp = fopen("step.info", "a+")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN step.info FILE. \n");
      exit(0);
    }
  
  fprintf(fp,"(%d) done photon %.2e out of %.2e --- %3.0f %s ---  TIME: %.2f min  dT= %.2f min\n",s,
	  photon_i, 
	  (1*npackets), 
	  100*photon_i/npackets, "%",
	  rtime0+(time(NULL)-time0)/60, 
	  (time(NULL)-(*time1_ptr))/60);
  *time1_ptr=time(NULL);

  fclose(fp);
	 
  return;
}

/********* final update step info file **********************************************/

void final_update_step_info(photon_i,npackets,time0, time1_ptr,s,rtime0)
     double photon_i;
     double npackets;
     double time0;
     double *time1_ptr;
     int s;
     double rtime0;
{
  FILE *fp;

  if (( fp = fopen("step.info", "a+")) == NULL)
   {
       fp = fopen("step.info", "w");
   }
  
  (*time1_ptr)=time(NULL);

  fprintf(fp,"(%d) done photon %.2e out of %.2e --- %3.0f %s ---  TIME: %.2f min\n",s, 
	  photon_i, (1*npackets), 100*photon_i/npackets,"%", rtime0+((*time1_ptr)-time0)/60);
  fclose(fp);
  
  return;
}


void write_lambda_opacity(fbin,opacity)
     struct fbin_type fbin;
     struct opacity_type opacity;
{
  int i;
  FILE *fp;

  if (( fp = fopen("opacity.extra.dat", "w")) == NULL)
    {fprintf(stderr, "FAILED TO OPEN opacity.extra FILE. \n");return;}
  
  for (i =0; i<=NFBINS; i++)   
    
    fprintf(fp,"%7.3f %.4e \n",
	    fbin.lambda[i], 
	    (k_abs(c_CONST/(fbin.lambda[i]*1e-4),opacity)+sigma_scat(c_CONST/(fbin.lambda[i]*1e-4),opacity)));
 
  fclose(fp);

  return;

} 




