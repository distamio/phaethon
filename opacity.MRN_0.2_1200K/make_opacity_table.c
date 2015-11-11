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


struct opacity_type opacity;

int main()

{
  int i;
  double junk;
  double lambda[OPADIM],k1[OPADIM],s3[OPADIM],k3[OPADIM];

FILE *fp;

   /* Open  file */
   if (( fp = fopen("opacity.dat_no_extra", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat_no_extra FILE. \n");
         return;
         }

for (i = 0;i <98; i++)
         
fscanf(fp,"%lf %lf %lf\n",

       &lambda[i],&k3[i],&s3[i]);
    
   fclose(fp);

for (i = 0;i <98; i++)
  {

    opacity.lambda[i]=lambda[i];
    opacity.abs[i]=k3[i];
    opacity.scat[i]=s3[i];
  }


/* FIND OPACITY >1000 microns..*/  

/* 1. continue as 1/lambda  */

for (i =98;i <OPADIM; i++)

  {
    opacity.lambda[i]=lambda[97]+(i-97)*(1.e6-lambda[97])/(OPADIM-97.);
    opacity.abs[i]=opacity.abs[97]*lambda[97]/opacity.lambda[i];
    opacity.scat[i]=0.;
  }


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
