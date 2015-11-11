/*
   Calculates Planck mean opacity in cgs units.

   compile: cc -o mopacity mopacity.c opacity.c ../main/inout.c ../main/bb.c -lm 


  [input:  opacity.dat (table with the opacities)]
   output: mopacity.dat (Temp vs planck mean opacity (cgs))

  opacity.dat file should contain lambda(microns), absorption, scattering    
                               (ascending order)

   NOTE: MAY NEEDS TO BE CHANGED FOR DIFFERENT OPACITY BECAUSE THE INTEGRATION
         LIMITS MAYBE DIFFERENT

*/ 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../main/structs.h"
#include "../opacity/params.opa.h"
#include "../params.h"
#include "../main/constants.h"
#include "../main/functions.h"
#include "../main/inout.h"
#include "../params.h"
#include "../main/opacity.h"

struct opacity_type opacity;
struct mean_op_type mean_op;

int main()

{
int i;
double f1,f2,a, a_aux,c2,h3,k2;
double lambda1,lambda2,aux;

/* read the opacity table (opacity should be read in cgs units !!!) */
/* opacity table should be lambda, abs, scat, increasing lambda     */

read_opacity(&opacity); 

/* calculate first and last frequency of the table */

f1=c_CONST/(1e-4*opacity.lambda[OPADIM-1]);
f2=c_CONST/(1e-4*opacity.lambda[0]);

/* some constants */

c2=pow(c_CONST,2.);
h3=pow(h_CONST,3.);
k2=pow(k_CONST,2);

a_aux=(h_CONST/k_CONST)*10*f2;

aux=1.*(MAX_TEMP-1)/(TEMPDIM-1);

for (i=1;i<TEMPDIM;i++)

  {

   mean_op.temp[i]=aux*i;

   printf("T= %.3f K\n",  mean_op.temp[i]);

   a=a_aux/mean_op.temp[i];
 
   /* i calculate the integral only were the BB function is important

   ie. 0<x<40  (maybe too much already) --> 0<freq<8.e11*Temp */

   mean_op.opacity[i]=M_PI*
     (simp_kb(f1/100.,8.e11*mean_op.temp[i],mean_op.temp[i]))   
   /(sigma_CONST*pow(mean_op.temp[i],4.));  
  }
  mean_op.opacity[0]=mean_op.opacity[1];

/* output mopacity.dat file (temp, mean opacity (cgs)) */

write_mean_opacity(mean_op);

return 0;
}


/******* Function calculate integral (k black) (using Simpson's rule) ********/
 
double simp_kb(a,b,temp) 
 
     double a; /* initial frequency */
     double b; /* final frequency   */
     double temp; 

{ 
double n=16; 
double sum=0.;
double space; 
int counter;   
  
space=(b-a)/pow(2.,n);
 
sum=black(a,temp)*k_abs(a,opacity);  

for (counter=1 ; counter<=pow(2.,n)-2 ; counter=counter+2.) 
        
sum=sum+4.*(black(a+space*counter,temp)*k_abs(a+space*counter,opacity))
      +2.*(black(a+space*(counter+1),temp)*k_abs(a+space*(counter+1),opacity));
  
sum=sum
  +4.*(black(a+space*(pow(2.,n)-1),temp)*k_abs(a+space*(pow(2.,n)-1),opacity))+
  +black(b,temp)*k_abs(b,opacity);
 
sum=sum*(space/3.); 

return(sum); 
} 

