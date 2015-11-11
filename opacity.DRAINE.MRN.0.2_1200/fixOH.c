/* Makes an opacity table from OH opacities

   with the same params.opa.dat as the Draine2003

  cc -g -o fixOH fixOH.c ../main/inout.c  -lm

  output: opacity2.dat (lambda (microns, abs, scat)

*/
#define MAX_LINE     50
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../opacity/params.opa.h"
#include "../main/constants.h"


#define DUST_TO_GAS  (1/124.)
#define OPADIM_IN 500 /* dim of the initial opa file */
#define OPADIM_DR 553

/* assumes OPADIM_DR >OPADIM_IN */

struct opacity_type
       {
        double lambda[OPADIM_IN];
        double abs[OPADIM_IN];
        double scat[OPADIM_IN];
        double albedo[OPADIM_IN];
       };

struct opacity_type opacity;

double k_abs(double,struct opacity_type);
double sigma_scat(double,struct opacity_type);
   

main()
{
  int i,k; 

  double lambda[OPADIM_DR],newlambda[OPADIM_DR],k1[OPADIM_DR],k2[OPADIM_DR],k3[OPADIM_DR];
  double total[OPADIM_DR],albedo[OPADIM_DR];
  char buff[MAX_LINE];
  
  FILE *fp;
  

/* read OH opacity file opacity.in.dat */
    
  if (( fp = fopen("opacity.in.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
      exit(0);
    }
  
  for (i = 0;i <OPADIM_IN; i++)
    {    
      fscanf(fp,"%lf %lf %lf\n",&opacity.lambda[i],&opacity.abs[i],& opacity.scat[i]);
      
      opacity.albedo[i]= opacity.scat[i]/(opacity.scat[i]+opacity.abs[i]);
    }  
  
  fclose(fp);


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
       
      i=i+1;
    }    

  fclose(fp);
 

 for (i =0; i<OPADIM_DR; i++)
   { 
     newlambda[i]=lambda[OPADIM_DR-1-i];    
   }
 

/* Open  file */
   if (( fp = fopen("opacity.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
         return;
	 }


for (i = 0;i<OPADIM_DR; i++)
     
  fprintf(fp,"%4e %4e %4e\n",newlambda[i],k_abs(c_CONST/(newlambda[i]*1e-4),opacity),
	  sigma_scat(c_CONST/(newlambda[i]*1e-4),opacity));
      
   fclose(fp);


return 0;

}






/********** function opacity k_abs *****************************************/

double k_abs(f,opacity)

 double f;
 struct opacity_type opacity;

{
 double lambda;
 int index;
 int k;

/* find photon lambda in microns */

 lambda=(c_CONST/f)*1e4;


 if (lambda<=opacity.lambda[0])  return(opacity.abs[0]);

 else
     
      if (lambda>=opacity.lambda[OPADIM_IN-1])  return(opacity.abs[OPADIM_IN-1]);
  
      else 

	{
         index=binary_search_opacity(lambda,0,OPADIM_IN-1,opacity);

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
 int k;



 

/* find photon lambda in microns */

 lambda=(c_CONST/f)*1e4;


 if (lambda<=opacity.lambda[0])  return(opacity.scat[0]);

 else
     
      if (lambda>=opacity.lambda[OPADIM_IN-1])  return(opacity.scat[OPADIM_IN-1]);
  
      else 

	{
         index=binary_search_opacity(lambda,0,OPADIM_IN-1,opacity);

         return( (opacity.scat[index+1]*(lambda-opacity.lambda[index])+ 
		  (opacity.scat[index]*(opacity.lambda[index+1]-lambda)) )/
		 (opacity.lambda[index+1]-opacity.lambda[index]));
	}

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

