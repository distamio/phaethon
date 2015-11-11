/* Makes an opacity table from OH opacities

  cc -g -o fixd fixdraine.c ../main/inout.c opacity.c -lm

  output: opacity2.dat (lambda (microns, abs, scat)

*/
#define MAX_LINE     50
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

#define DUST_TO_GAS  (1/124.)

#define OPADIM_DR  OPADIM

struct opacity_type opacity;
struct opacity_type opacity2;


main()
{
  int i,k;

  double lambda[OPADIM_DR],k1[OPADIM_DR],k2[OPADIM_DR],k3[OPADIM_DR];
  double total[OPADIM_DR],albedo[OPADIM_DR];
  char buff[MAX_LINE];
  
  FILE *fp;
  
  /* Open  file */
  if (( fp = fopen("draine03.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN draine03.dat FILE. \n");
      return;
    }
  
  i=0;

  for (k = 0; k<800 & getdata(fp, buff); k++)
    {         
     
      sscanf(buff,"%lf %lf %lf %lf %lf %lf \n",
	     &lambda[i],&albedo[i],&k1[i],&k2[i],&total[i],&k3[i]);
       
      /* printf("%lf %d %lf \n", lambda[i], i, total[i]); */ 
      i=i+1;
    }    

  fclose(fp);
  
 
  
  for (i =0; i<OPADIM_DR; i++)
    { 
      opacity.lambda[i]=lambda[OPADIM_DR-1-i];
      opacity.abs[i]=(1-albedo[OPADIM_DR-i-1])*total[OPADIM_DR-i-1];
      opacity.scat[i]=albedo[OPADIM_DR-i-1]*total[OPADIM_DR-i-1];
      printf("%lf %d %lf %lf %lf\n ", lambda[OPADIM_DR-i-1], i, opacity.lambda[i], opacity.abs[i], opacity.scat[i]);
    }
  
  



/* Open  file */
   if (( fp = fopen("opacity2.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity2.dat FILE. \n");
         return;
	 }


for (i = 0;i <OPADIM_DR; i++)
     
  fprintf(fp,"%4e %4e %4e\n",opacity.lambda[i],opacity.abs[i]*DUST_TO_GAS,opacity.scat[i]*DUST_TO_GAS);
      
   fclose(fp);


return 0;

}







