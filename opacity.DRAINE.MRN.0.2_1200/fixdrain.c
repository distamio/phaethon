/* Makes an opacity table from OH opacities

  cc -g -o fixd fixdraine.c opacity.c -lm

  output: opacity2.dat (lambda (microns, abs, scat)

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

#define DUST_TO_GAS  124

#define OPADIM_DR  552

struct opacity_type opacity;
struct opacity_type opacity2;

double sigma_scat_special(double,struct opacity_type);

main()
{
  int i;
  double density,junk[30];
  double lambda[OPADIM_DR],k1[OPADIM_DR],k2[OPADIM_DR],k3[OPADIM_DR];
  double total[OPADIM_DR];
 
FILE *fp;

   /* Open  file */
   if (( fp = fopen("draine03.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN draine03.dat FILE. \n");
         return;
         }


for (i = 0;i <600 && getdata(fp, buff); i++)
         
fscanf(fp,"%lf %lf %lf %lf %lf %lf \n",

       &lambda[i],&albedo[i],&k1[i],&k2[i],&total[i],&k3[i]);
    
   fclose(fp);

for (i =0; i<OPADIM_DR; i--)
  {
    opacity.lambda[i]=lambda[OPADIM_DR-i];
    opacity.abs[i]=(1-albedo[OPADIM_DR-i])*total[OPADIM_DR-i];
    opacity.scat[i]=albedo[OPADIM_DR-i]*total[OPADIM_DR-i];
  }





/* Open  file */
   if (( fp = fopen("opacity2.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity2.dat FILE. \n");
         return;
	 }


for (i = 0;i <OPADIM_DR; i++)
     
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
