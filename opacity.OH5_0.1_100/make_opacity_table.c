/* Makes an opacity table from OH opacities

  cc -g -o make_opacity_table make_opacity_table.c opacity.c -lm

  output: opacity.dat (lambda (microns, abs, scat)---- scat=MRN_scat/2.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "params.opa.h"
#include "constants.h"
#include "functions.h"
#include "inout.h"
#include "opacity.h"


#define OPADIM_MRN  552 
#define OPADIM_OH    66
#define OPADIM_TEMP 500

#define GAS_TO_DUST 100
#define SCATTER_FACTOR 0.5 
struct opacity_type opacity;
struct opacity_type opacity2;
struct opacity_type opacity3;

double sigma_scat_special(double,struct opacity_type,int);
double sigma_abs_special(double,struct opacity_type,int);
main()
{
  int i;
  double density,junk[30];
  double lambda[OPADIM_TEMP],k1[OPADIM_TEMP],k2[OPADIM_TEMP],k3[OPADIM_TEMP];
  double lambda_MRN[OPADIM_MRN],k_MRN[OPADIM_MRN],s_MRN[OPADIM_MRN];
  double dlambda;
  
FILE *fp;




int k;


// make lambda table for opacity (in microns)

opacity3.lambda[0]=L1; /* starting lambda (microns) */
dlambda=log10(L2/L1)/OPADIM; /* logarithic fbins      */

for (i=1;i<=OPADIM;i++)
  {
   opacity3.lambda[i]=pow(10,(log10(opacity3.lambda[i-1])+dlambda));
  }

// create temporary opacity table depending on the give table

   /* Open  file */
   if (( fp = fopen("opacity.dat_OH5", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_OH5 FILE. \n");
         return;
         }

for (i = 0;i <OPADIM_OH; i++)
         
fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf \n",

       &density,&lambda[i],&k1[i],&lambda[i],&k2[i],&lambda[i],&k3[i]);
    
   fclose(fp);

for (i = 0;i <OPADIM_OH; i++)
  {
    opacity.lambda[i]=lambda[i];
    opacity.abs[i]=k3[i];
    opacity.scat[i]=0.;
  }


/*  continue OH5  to longer wavelengthts*/

for (i =OPADIM_OH;i <OPADIM_TEMP-129; i++)

   {
     opacity.lambda[i]=lambda[OPADIM_OH-1]+(i-(OPADIM_OH-1))*(1.e5-lambda[OPADIM_OH-1])/(OPADIM_TEMP-OPADIM_OH);
     opacity.abs[i]=pow(10,5.497-1.78*log10(opacity.lambda[i]))/k3[OPADIM_OH-1];
     opacity.scat[i]=0.;
   }

/* move everything 129 positions up in order to put MRN  and arrange

   for the dust to gas ratio */

for (i =OPADIM_TEMP-130; i>=0; i--)

   {
     opacity.lambda[i+129]= opacity.lambda[i];
     opacity.abs[i+129]=opacity.abs[i]/GAS_TO_DUST;
     opacity.scat[i+129]= opacity.scat[i];
   }

/* extend to 0.1 microns with MRN ...*/

/* Open  MRN file */

   if (( fp = fopen("opacity.dat_MRN_draine03", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_MRN_draine03 FILE. \n");
         return;
         }

 for (i = 0;i <OPADIM_MRN ; i++)

 {
 fscanf(fp,"%lf %lf %lf  \n",

        &lambda_MRN[i],&k_MRN[i],&s_MRN[i]);
 
 }
 fclose(fp);

 /* add MRN before OH5 */

 for (i = 0;i <129; i++) 

   {
     opacity.lambda[i]=lambda_MRN[i];
     opacity.abs[i]=k_MRN[i]/(k_MRN[129]/opacity.abs[129]);
   }

 /* calculate scattering  */


 /* Open  file */
if (( fp = fopen("opacity.dat_MRN_draine03", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_MRN_draine03 FILE. \n");
         return;
         }

 for (i = 0;i <OPADIM_MRN; i++)
 fscanf(fp,"%lf %lf %lf\n",&opacity2.lambda[i],&opacity2.abs[i],&opacity2.scat[i]);
      
 fclose(fp);

 
 for (i = 129 ;i <OPADIM_TEMP; i++)
 { 
  opacity.scat[i]=sigma_scat_special(opacity.lambda[i],opacity2,OPADIM_MRN);

 }

 for  (i =0;i<129; i++) 
 {
  opacity.scat[i]=s_MRN[i];

 }








// make new table

for (i = 0 ;i <OPADIM; i++)
 { 
  opacity3.scat[i]=sigma_scat_special(opacity3.lambda[i],opacity,OPADIM_TEMP);
  opacity3.abs[i]=sigma_abs_special(opacity3.lambda[i],opacity,OPADIM_TEMP);
 }




if (( fp = fopen("opacity.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
         return;
         }

for (i = 0; i <OPADIM; i++)
    fprintf(fp,"%.5e %.4e %.4e\n",opacity3.lambda[i],opacity3.abs[i],opacity3.scat[i]*SCATTER_FACTOR);
      
fclose(fp);


return 0;

}






/********** function opacity sigma_scat *************************************/

double sigma_scat_special(lambda,opacity,dim) 

double lambda;/* in microns */
struct opacity_type opacity;
int dim;

{
  int index;

 if (lambda<=opacity.lambda[0])  return(opacity.scat[0]);
 
 else

      if (lambda>=opacity.lambda[dim-1])  return(opacity.scat[dim-1]);

      else

        {
         index=binary_search_opacity(lambda,0,dim-1,opacity);

         return( (opacity.scat[index+1]*(lambda-opacity.lambda[index])+  
                  (opacity.scat[index]*(opacity.lambda[index+1]-lambda)) )/
                 (opacity.lambda[index+1]-opacity.lambda[index]));
        }

}



/********** function opacity sigma_abs *************************************/

double sigma_abs_special(lambda,opacity,dim)

double lambda;/* in microns */
struct opacity_type opacity;
int dim;

{
  int index;

 if (lambda<=opacity.lambda[0])  return(opacity.abs[0]);
 
 else

      if (lambda>=opacity.lambda[dim-1])  return(opacity.abs[dim-1]);

      else

        {
         index=binary_search_opacity(lambda,0,dim-1,opacity);

         return( (opacity.abs[index+1]*(lambda-opacity.lambda[index])+  
                  (opacity.abs[index]*(opacity.lambda[index+1]-lambda)) )/
                 (opacity.lambda[index+1]-opacity.lambda[index]));
        }

}

