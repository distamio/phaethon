
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "constants.h"
#include "params.h"
#include "phaethon.params.h"
#include "params.grid.h"
#include "grid.functions.h"

/***** density function *************************************************/

double dens(r,theta,phi) 
       
     double  r,theta,phi; 
{
  double c;

 if (DISC==1 && r<DISC_R) return(dens_disc(r,theta,phi));
 
 else
 
  { 
   if (r<R_DUST) return(0);

   else

   {   
    if (JET==1 && fabs(r*cos(theta))>fabs(JET_ALPHA*pow(r*sin(theta),JET_BETA))) return(JET_DENS);
    else
    {
    if (GEOMETRY==0)
     {
      if (r>R_MAX_CORE) return(DENS_EXT_FACTOR*dens(R_MAX_CORE,0.,0.));	    
      else  return(DENS0/(1+pow(r/R0,2.))); 
     }

    if (GEOMETRY==1)
     {
      c=r/R0;

      if (r>R_MAX_CORE) return(DENS_EXT_FACTOR*dens(R_MAX_CORE,0.,0.));	    
      else  return(DENS0*(1.+(A/49.)*c*c*pow(sin(theta),SIN_POWER))/((1+c*c)*(1+c*c)));
     }

    if (GEOMETRY>1)
     {
      printf("\n:: GEOMETRY>=2: not ready yet\n"); exit(0);
	
     }

     }
   }
 } 
}

/************* mass function *******************************************/

double mass(r1,r2,theta1,theta2,phi1,phi2) 

     double r1,r2,theta1,theta2,phi1,phi2; 
     
{ 
  double I1,I2,I3;

  if (DISC==1 && r2<=DISC_R)
 
   return(mass_disc(r1,r2,theta1,theta2,phi1,phi2));

  else

 {

if (r2<1.00000000000001*R_DUST) return(0.);

  else
{
    if (JET==1 && fabs(r2*cos(theta2))>fabs(JET_ALPHA*pow(r2*sin(theta2),JET_BETA)))  return(2*M_PI*JET_DENS*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);
    else

  {

  if (GEOMETRY==0) 
  {
  if (r1<R_MAX_CORE) return(DENS0*(phi2-phi1)*(cos(theta1)-cos(theta2))*R0*R0*R0*(((r2-r1)/R0)-atan(r2/R0)+atan(r1/R0)));
  else  return(DENS_EXT_FACTOR* dens(R_MAX_CORE,0.,0.)*(phi2-phi1)*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);
  }
  
 if (GEOMETRY==1) 
  {

   if (SIN_POWER==1) I3=((0.25*(sin(2*theta1)-sin(2*theta2)))+0.5*(theta2-theta1));
	      
	      else 
		{
		  if  (SIN_POWER==2) I3=0.75*(cos(theta1)-cos(theta2))-(1/12.)*(cos(3*theta1)-cos(3*theta2));
		  
		  else 
		    
		    if  (SIN_POWER==4)
		      I3=((5./8.)*(cos(theta1)-cos(theta2))-(5/48.)*(cos(3*theta1)-cos(3*theta2))   
			  +(1/80.)*(cos(5*theta1)-cos(5*theta2))); 
		}
	      
	      I1 =R0*R0*(-((-r1+R0*atan(r1/R0)+r1*(r1/R0)*atan(r1/R0))/(2*(1+(r1/R0)*(r1/R0))))+
			 ( (-r2+R0*atan(r2/R0)+r2*(r2/R0)*atan(r2/R0))/(2*(1+(r2/R0)*(r2/R0)))));
	      
	      I2=-((3*R0*R0*r1+2*r1*r1*r1-3*R0*R0*R0*atan(r1/R0)-3*R0*r1*r1*atan(r1/R0))/(2*(1+(r1/R0)*(r1/R0))))+
		((3*R0*R0*r2+2*r2*r2*r2-3*R0*R0*R0*atan(r2/R0)-3*R0*r2*r2*atan(r2/R0))/(2*(1+(r2/R0)*(r2/R0))));

  if (r1<R_MAX_CORE) return(DENS0*(phi2-phi1)*((cos(theta1)-cos(theta2))*I1+(1.*A/49.)*I3*I2));
  else  return(DENS_EXT_FACTOR*dens(R_MAX_CORE,0.,0.)*(phi2-phi1)*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);
  }
}
 }	 
} 
}
/***** density function *************************************************/

double dens_disc(r,theta,phi)

     double  r,theta,phi;
{
  double theta_opening; double theta0;

  if (r<0.9999999*R_DUST) return (0);
  else
    {
      
     if (JET==1 && fabs(r*cos(theta))>fabs(JET_ALPHA*pow(r*sin(theta),JET_BETA))) return(JET_DENS);

      else
        {
              if (theta<0.5*M_PI-DISC_THETA_OPENING*DISC_ALPHA || theta>0.5*M_PI+DISC_THETA_OPENING*DISC_ALPHA) 

              return(DISC_AMBIENT_DENS);
 	      else
{
          theta=0.5*M_PI-theta;
 		if (DISC_ALPHA==0) return (DISC_DENS0*pow((1*AU)/r,2.));
		else  
                return (DISC_DENS0*pow((1*AU)/r,2.)*exp(-0.5*pow(DISC_R/DISC_Z,2.)*theta*theta));
}         
        }
    }
}

/************* mass function *******************************************/

double mass_disc(r1,r2,theta1,theta2,phi1,phi2)  /* normal polar, azimuthal angles */

     double r1,r2,theta1,theta2,phi1,phi2;

{
  double alpha=DISC_ALPHA; double theta_opening; /* in degrees */ double theta0;


  if (fabs(theta1-theta2)<1e-15) return (0.);
  else
    {
  
      if (r2<=R_DUST)  return(0.*2*M_PI*DENS0*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);
      else
        {
          if (JET==1 && fabs(r2*cos(theta2))>fabs(JET_ALPHA*pow(r2*sin(theta2),JET_BETA)))  
return(2*M_PI*JET_DENS*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);

          else
            {

              if (theta2<0.5*M_PI-DISC_THETA_OPENING*DISC_ALPHA || theta1>0.5*M_PI+DISC_THETA_OPENING*DISC_ALPHA) 

		return(2*M_PI*DISC_AMBIENT_DENS*(r2*r2*r2-r1*r1*r1)*(cos(theta1)-cos(theta2))/3.);

              else
{
              theta1=0.5*M_PI-theta1;
              theta2=0.5*M_PI-theta2;

              if (DISC_ALPHA==0) 	
 	      return ((phi2-phi1)*pow(1*AU,2)*DISC_DENS0*(r2-r1));
     		else
              {
	       if (r1>DISC_INNER_GAP) 
              return ((phi2-phi1)*pow(1*AU,2)*DISC_DENS0*
                      (r2-r1)*alpha*1.253*(erf(0.707*theta1/alpha)-erf(0.707*theta2/alpha))
                      );
	       else 
              return ((phi2-phi1)*pow(1*AU,2)*DISC_DENS0*
                      (r2-DISC_INNER_GAP)*alpha*1.253*(erf(0.707*theta1/alpha)-erf(0.707*theta2/alpha))
                      );

	       }
            }
}
        }
    }
}

