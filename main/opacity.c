/*

   mopacity works for every T (not just integer T) [14.01.02]

   This file containts the functions that return the value of opacities for
   every frequency.

   the input in these functions is frequency f and opacity (just in case you 
   need to use a table with opacities)

   NOTE: THIS SHOULD BE CHANGED FOR DIFFERENT OPACITIES !!!

   D. Stamatellos [19.11.01] 

*/

#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "structs.h"
#include "opacity.h"



/********** function opacity k_abs *****************************************/

double k_abs(f,opacity)

 double f;
 struct opacity_type opacity;

{
 double lambda;
 int index;

/* find photon lambda in microns */

 lambda=(c_CONST/f)*1e4;


 if (lambda<=opacity.lambda[0])  return(opacity.abs[0]);

 else
     
      if (lambda>=opacity.lambda[OPADIM-1])  return(opacity.abs[OPADIM-1]);
  
      else 

	{
         index=binary_search_opacity(lambda,0,OPADIM-1,opacity);

         return( (opacity.abs[index+1]*(lambda-opacity.lambda[index])+ 
		  (opacity.abs[index]*(opacity.lambda[index+1]-lambda)) )/
		 (opacity.lambda[index+1]-opacity.lambda[index]));
	}

}



/********** function opacity sigma_scat *************************************/

double sigma_scat(f,opacity)

double f;
struct opacity_type opacity;

{
 double lambda;
 int index;

/* find photon lambda in microns */

 lambda=(c_CONST/f)*1e4;


 if (lambda<=opacity.lambda[0])  return(opacity.scat[0]);

 else
     
      if (lambda>=opacity.lambda[OPADIM-1])  return(opacity.scat[OPADIM-1]);
  
      else 

	{
         index=binary_search_opacity(lambda,0,OPADIM-1,opacity);

         return( (opacity.scat[index+1]*(lambda-opacity.lambda[index])+ 
		  (opacity.scat[index]*(opacity.lambda[index+1]-lambda)) )/
		 (opacity.lambda[index+1]-opacity.lambda[index]));
	}

}


/********** function albedo *****************************************/

double albedo(f,opacity)

double f;
struct opacity_type opacity;

{
 return(sigma_scat(f,opacity)/(k_abs(f,opacity)+sigma_scat(f,opacity)));
}

/**************** function mean_opacity  ******************************/

double mean_opacity(temp,mean_op_ptr)

     double temp; 
     struct mean_op_type * mean_op_ptr;
{
 
  int index;
 int middle,low,high;

 low=0;
 high=TEMPDIM-1;

 while (low<=high) {

 middle= (low+high)/2;

 if (temp>=(*mean_op_ptr).temp[middle] && temp<(*mean_op_ptr).temp[middle+1]) 

 
   break;
 
 else if (temp<(*mean_op_ptr).temp[middle])
   
   high=middle-1;

 else low=middle+1;
 
}

index=middle;

return( ((*mean_op_ptr).opacity[index+1]*(temp-(*mean_op_ptr).temp[index])+ 
		  ((*mean_op_ptr).opacity[index]*((*mean_op_ptr).temp[index+1]-temp)) )/
		 ((*mean_op_ptr).temp[index+1]-(*mean_op_ptr).temp[index]));

}



/******************** modified binary search **************************/

int binary_search_opacity (lambda, low, high, opacity)


double lambda;
int low;
int high;
struct opacity_type opacity;

{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (lambda>=opacity.lambda[middle] && lambda<opacity.lambda[middle+1]) 

    return middle;

 else if (lambda<opacity.lambda[middle])

         high=middle-1;

 else low=middle+1;

 }

return -1; /* search key not found */
}



/***************** function find_new_f *************************************/

/* finds the frequency bin  of the reemitted photon packet */

void find_photon_fpoint(p,i,reemit_ptr)
  
 struct photon_packet_type p[];
 int  i;
 struct reemit_type * reemit_ptr;

{
 int low=0;
 int  high=NFREQ-1;
 int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (p[i].f>=(*reemit_ptr).freq[middle] && p[i].f<(*reemit_ptr).freq[middle+1]) 

    break;

 else if (p[i].f<(*reemit_ptr).freq[middle])

         high=middle-1;

 else low=middle+1;

 }

 p[i].fpoint=middle;

return;

}

/********* create a table that has the opacities for each fpoint ********************/

void  create_opacity_table(opa_ptr ,reemit_ptr ,opacity)
     
     struct opacity_table * opa_ptr;
     struct reemit_type  *reemit_ptr;
     struct opacity_type opacity;
{
  int i;
  
  for (i = 0; i <NFREQ; i++)
    {     
      (*opa_ptr).abs[i]=k_abs((*reemit_ptr).freq[i],opacity);
      (*opa_ptr).scat[i]=sigma_scat((*reemit_ptr).freq[i],opacity);
      (*opa_ptr).albedo[i]= (*opa_ptr).scat[i]/((*opa_ptr).abs[i]+(*opa_ptr).scat[i]);
    }
  return;
}
  
  
