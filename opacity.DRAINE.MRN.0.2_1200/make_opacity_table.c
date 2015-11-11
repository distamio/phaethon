/* Makes an opacity table from OH opacities

  cc -g -o make_opacity_table make_opacity_table.c opacity.c -lm

  output: opacity.dat (lambda (microns, abs, scat)/ relative to 1micron absorptive opacity

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../main/structs.h"
#include "../opacity/params.opa.h"
#include "../main/constants.h"
#include "../main/functions.h"
#include "../main/inout.h"
#include "../main/opacity.h"

#define DUST_TO_GAS  100

#define OPADIM_MRN   97

struct opacity_type opacity;
struct opacity_type opacity2;

double sigma_scat_special(double,struct opacity_type);

main()
{
  int i;
  double density,junk[30];
  double lambda[OPADIM],k1[OPADIM],k2[OPADIM],k3[OPADIM];
  double lambda_MRN[OPADIM],k_MRN[OPADIM],s_MRN[OPADIM];
 
FILE *fp;

   /* Open  file */
   if (( fp = fopen("opacity.dat_OH", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_OH FILE. \n");
         return;
         }

for (i = 0;i <67; i++)
         
fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf \n",

       &density,&lambda[i],&k1[i],&lambda[i],&k2[i],&lambda[i],&k3[i]);
    
   fclose(fp);

for (i = 0;i <67; i++)
  {
    opacity.lambda[i]=lambda[i];
    opacity.abs[i]=k3[i];
    opacity.scat[i]=0.;
  }

/* FIND OPACITY >1000 microns.. 2 ways---> */  

/* 1. continue OH5 as 1/lambda  */
/*
for (i =67;i <OPADIM; i++)

  {
    opacity.lambda[i]=lambda[66]+(i-66)*(1.e6-lambda[66])/(OPADIM-67.);
    opacity.abs[i]=opacity.abs[66]*lambda[66]/opacity.lambda[i];
    opacity.scat[i]=0.;
  }

*/

/* 2. continue OH5 (need to modify...opacity) */

for (i =67;i <OPADIM-19; i++)

   {
     opacity.lambda[i]=lambda[66]+(i-66)*(1.e6-lambda[66])/(OPADIM-67.);
     opacity.abs[i]=pow(10,5.497-1.78*log10(opacity.lambda[i]))/k3[66];
     opacity.scat[i]=0.;
   }

/* move everything 19 positions up in order to put MRN  and arrange

   for the dust to gas ratio */

for (i =OPADIM-20; i>=0; i--)

   {
     opacity.lambda[i+19]= opacity.lambda[i];
     opacity.abs[i+19]=opacity.abs[i]/DUST_TO_GAS;
     opacity.scat[i+19]= opacity.scat[i];
   }

/* extend to 0.1 microns with MRN ...*/

/* Open  MRN file */

   if (( fp = fopen("opacity.dat_MRN", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_MRN.dat FILE. \n");
         return;
         }

 for (i = 0;i <97; i++)

 fscanf(fp,"%lf %lf %lf  \n",

        &lambda_MRN[i],&k_MRN[i],&s_MRN[i]);
     
 fclose(fp);

 /* add MRN before OH5 */

 for (i = 0;i <19; i++) 

   {
     opacity.lambda[i]=lambda_MRN[i];
     opacity.abs[i]=k_MRN[i]/(k_MRN[18]/opacity.abs[19]);
   }


 /* calculate scattering  */


 /* Open  file */
   if (( fp = fopen("opacity.dat_MRN", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_MRN FILE. \n");
         return;
         }

  for (i = 0;i <OPADIM_MRN; i++) 

  fscanf(fp,"%lf %lf %lf\n",&opacity2.lambda[i],&opacity2.abs[i],&opacity2.scat[i]);
      
  fclose(fp);

 for (i = 0;i <97; i++) 

  opacity.scat[i]=sigma_scat_special(opacity.lambda[i],opacity2)/(k_MRN[18]/opacity.abs[19]);

 for  (i = 97;i<OPADIM; i++) opacity.scat[i]=0;
  

/* Open  file */
   if (( fp = fopen("opacity.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
         return;
	 }


for (i = 0;i <OPADIM; i++)
     
      fprintf(fp,"%.5e %4e %4e\n",opacity.lambda[i],opacity.abs[i],opacity.scat[i]);
      
   fclose(fp);


return 0;

}








/********** function opacity sigma_scat *************************************/

double sigma_scat_special(lambda,opacity)

double lambda;/* in microns */
struct opacity_type opacity;

{
  int index;

 if (lambda<=opacity.lambda[0])  return(opacity.scat[0]);

 else
     
      if (lambda>=opacity.lambda[OPADIM_MRN-1])  return(opacity.scat[OPADIM_MRN-1]);
  
      else 

	{
         index=binary_search_opacity(lambda,0,OPADIM_MRN-1,opacity);

         return( (opacity.scat[index+1]*(lambda-opacity.lambda[index])+ 
		  (opacity.scat[index]*(opacity.lambda[index+1]-lambda)) )/
		 (opacity.lambda[index+1]-opacity.lambda[index]));
	}

}
