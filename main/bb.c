/* This file contains blackbody functions and routines for calculating related
   integrals
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "structs.h"
#include "functions.h"
#include "constants.h"




/******** Function calculate integral_nblack *****************************/
 
/* finds the integral from a to b of the normalized bb function 

   (using Simpson's rule)    */

/* bb_const=sigma * T**4 /pi */

double simp_b(a, b, bb_const) 
 
 double a;
 double b;
 double bb_const; /* const to use if I want to convert to cgs */
{
  double n=10; 
  double sum=0.;
  double space; 
  int counter;  
  
  space=(b-a)/pow(2.,n);
 
  sum=nblack(a,bb_const);  

  for (counter=1 ; counter<=pow(2.,n)-2 ; counter=counter+2.) 
        
     sum=sum +4.*nblack(a+space*counter,bb_const)
             +2.*nblack(a+space*(counter+1.),bb_const);
  

  sum=sum+4.*nblack(a+space*(pow(2.,n)-1),bb_const)+nblack(b,bb_const);
 
  sum=sum*(space/3.); 

  return(sum); 
} 

/******** Normalized BlackBody Function *******************************/

double nblack(x,b_const) 

       double x;
       double b_const;

{
if (x<0.00003) 

   return ( b_const*(15./pow(M_PI,4))*pow(x,2));

else  

   return( b_const*(15./pow(M_PI,4))*(pow(x,3)/(exp(x)-1)) );
}

/******** BlackBody Function *******************************/

/* cgs units, frequency input*/

double black(f,t) 

     double f; /* frequency */
     double t; /* temperature */

{
 double x;
 double constant;

 x=(h_CONST/k_CONST)*(f/t);

 constant=2*h_CONST*pow(f,3)/(c_CONST*c_CONST);

if (x<0.00003) 

   return (constant/x);

else  

   return(constant/(exp(x)-1));

}

/******** BlackBody Function *******************************/

/* cgs units, wavelenth input */

double lblack(l,t) /**** this is B_lambda  **********/
  
     double l; /* wavelength (cm) */
     double t; /* temperature     */
     
{
  double x;
  double constant;
  
  x=(h_CONST/k_CONST)*(c_CONST/(l*t));
  
  constant=2*h_CONST*(c_CONST*c_CONST)/pow(l,5);
  
  if (x<=0.00003)  return (constant/x);
  
  if (x>=50 && x<100) return(constant*exp(-x));
  
  if (x>=100) return(0.);
  
  if (x>0.00003 && x<50)	return(constant/(exp(x)-1));
 
}




/******** BlackBody Function *******************************/

/* cgs units, wavelenth input */

double lblack_nu(l,t) /**** this is B_nu *********/

     double l; /* wavelength (cm) */
     double t; /* temperature     */

{
 double x;
 double constant;

x=(h_CONST/k_CONST)*(c_CONST/(l*t));

constant=2*h_CONST*(c_CONST/pow(l,3));

 if (x<=0.00003)  return (constant/x);

 if (x>=50 && x<100) return(constant*exp(-x));

 if (x>=100) return(0.);

 if (x>0.00003 && x<50)	return(constant/(exp(x)-1));

}


/******** d(BlackBody Function)/dT *******************************/

/* cgs units */

double dblack(f,t) 

       double f;
       double t;

{
 double x;
 double constant1,constant2;

 x=(h_CONST/k_CONST)*(f/t);

 constant1=2*((h_CONST*h_CONST)/(k_CONST*c_CONST*c_CONST));
 constant2=2*k_CONST/(c_CONST*c_CONST);

if (x<0.00003) 

   return (constant2*pow(f,2));

if (x>100) 

   return (0);
 
 if (x<100 && x>0.00003)   

   return(constant1*pow(f,4)*pow(t,-2.)*(exp(x)/pow((exp(x)-1),2)));

}

/******************** modified binary search **************************/

int binary_search_isrf (lambda, low, high, bisrf_ptr)


double lambda;
int low;
int high;
struct btype *bisrf_ptr;

{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (lambda>=(*bisrf_ptr).lambda[middle] && lambda<(*bisrf_ptr).lambda[middle+1]) 

    return middle;

 else if (lambda<(*bisrf_ptr).lambda[middle])

         high=middle-1;

 else low=middle+1;

 }

return -1; /* search key not found */
}

/********** function isrf *****************************************/

double isrf(lambda,bisrf_ptr)  /* returns Inu */

 double lambda; /* microns */
 struct btype *bisrf_ptr;

{
 int index;
 int isrf_dim=(*bisrf_ptr).dim;

 if (lambda<=(*bisrf_ptr).lambda[0])  return(0);

 else
     
      if (lambda>=(*bisrf_ptr).lambda[isrf_dim-1])  return((*bisrf_ptr).intensity[isrf_dim-1]);
  
      else 

	{
         index=binary_search_isrf(lambda,0,isrf_dim-1,bisrf_ptr);

        return( ((*bisrf_ptr).intensity[index+1]*(lambda-(*bisrf_ptr).lambda[index])+ 
		  ((*bisrf_ptr).intensity[index]*((*bisrf_ptr).lambda[index+1]-lambda)) )/
		 ((*bisrf_ptr).lambda[index+1]-(*bisrf_ptr).lambda[index]));


	}

}





/******** Function calculate integral_isrf *****************************/
 
/* finds the integral from a to b of the Black (1994) ISRF function or any other function eg. andre

   (using Simpson's rule)    */

double simp_isrf(a, b, bisrf_ptr,n) 
 
     double a; /* initial frequency to start integration */
     double b; /* final frequency  ... */
     struct btype *bisrf_ptr; 
     int n;
{
  double sum=0.;
  double space; 
  int counter;  
 

  space=(b-a)/pow(2.,n);
 
  sum=isrf(1.e4*c_CONST/a,bisrf_ptr);  

  for (counter=1 ; counter<=pow(2.,n)-2 ; counter=counter+2.) 
        
     sum=sum +4.*isrf(1.e4*c_CONST/(a+space*counter),bisrf_ptr)
             +2.*isrf(1.e4*c_CONST/(a+space*(counter+1.)),bisrf_ptr);
  

  sum=sum+4.*isrf(1.e4*c_CONST/(a+space*(pow(2.,n)-1)),bisrf_ptr)+isrf(1.e4*c_CONST/b,bisrf_ptr);
 
  sum=sum*(space/3.); 

  return(sum); 
} 
/******** Function calculate integral_star_sed *****************************/
 
/* finds the integral from a to b of the Black (1994) ISRF function or any other function eg. andre

   (using Simpson's rule)    */

double simp_star_sed(a, b, star,s,n) 
 
     double a; /* initial frequency to start integration */
     double b; /* final frequency  ... */
     struct star_type star[]; 
     int s;
     int n;
{
  double sum=0.;
  double space; 
  int counter;  
 
  space=(b-a)/pow(2.,n);
 
  sum=star_sed(1.e4*c_CONST/a,star,s);  

  for (counter=1 ; counter<=pow(2.,n)-2 ; counter=counter+2.) 
      
    {  
      sum=sum +4.*star_sed(1.e4*c_CONST/(a+space*counter),star,s)
	+2.*star_sed(1.e4*c_CONST/(a+space*(counter+1.)),star,s);
    
    }
  sum=sum+4.*star_sed(1.e4*c_CONST/(a+space*(pow(2.,n)-1)),star,s)+star_sed(1.e4*c_CONST/b,star,s);
 
  sum=sum*(space/3.); 

  return(sum); 
} 




/********** function star_sed*****************************************/

double star_sed(lambda,star,s)  /* returns Inu */

 double lambda; /* microns */
 struct star_type star[];
int s;

{
 int index;
 int isrf_dim=star[s].sed_dim;

 if (lambda<=star[s].lambda[0])  return(0);
 else
   if (lambda>=(star[s].lambda[isrf_dim-1]))  return(star[s].intensity[isrf_dim-1]);
      else 
	{
	  index=binary_search_star_sed(lambda,0,isrf_dim-2,star,s);

         return( (star[s].intensity[index+1]*(lambda-star[s].lambda[index])+ 
		  (star[s].intensity[index]*(star[s].lambda[index+1]-lambda)) )/
		 (star[s].lambda[index+1]-star[s].lambda[index]));
	}
}

/******************** modified binary search **************************/

  int binary_search_star_sed (lambda, low, high, star,s)


  double lambda;
 int low;
 int high;
 struct star_type star[];
 int s;
 
{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (lambda>=star[s].lambda[middle] && lambda<star[s].lambda[middle+1]) 

    return middle;

 else if (lambda<star[s].lambda[middle])

         high=middle-1;

 else low=middle+1;

 }

return -1; /* search key not found */
}
