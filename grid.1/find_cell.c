#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "structs.h"
#include "params.h"
#include "phaethon.params.h"
#include "params.grid.h"
#include "grid.functions.h"
#include "constants.h"
#include "phaethon.functions.h"

int binary_search_cell (double, int, int,ctypeptr);

/***** function find_photon_cell **************************************/

/* returns the index of the cell that the photon packet is in. If the 
   photon packet is outside the grid returns -1 */



int find_photon_cell(p,i,c,star,star_id)

 struct photon_packet_type p[];
 long i;
 struct cell_type c[];
 struct star_type star[];
 int star_id;

{
 int s,k,j,flag=0;
 int index=-100;
 double theta,aux,phi,theta_max;
 int flag1=0,flag2=0;
 double alpha,beta,mfp;
 
// printf("\nstar_id].z=%.4e c[NRCELLS-1].x=%.4e %.4e %.4e %.4e\n  ",star[star_id].z/AU,c[NRCELLS-1].x/AU,
// asin(c[NRCELLS-1].x/star[star_id].z)*180/M_PI,
// M_PI-asin(c[NRCELLS-1].x/star[star_id].z)*180/M_PI,
// acos(p[i].kz)*180/M_PI);
// exit(0);	 
   

 if ((star[star_id].type==3 || star[star_id].type==5) && p[i].flag==0 && star[star_id].R>R_MAX)
   {
	   
	   theta_max=asin(c[NRCELLS-1].x/star[star_id].R); // theta_max in (0,pi/2)          
	   theta=acos((-p[i].x*p[i].kx-p[i].y*p[i].ky-p[i].z*p[i].kz)/fabs(p[i].r)); // theta in (0,pi)
//	   printf("theta: %.f theta_max: %.f\n",theta*180/M_PI,theta_max*180/M_PI);
	   
     if (theta>theta_max) index=-1; 
     else
       { 
		
		mfp=fabs(p[i].r)*cos(theta)-pow(c[NRCELLS-1].x*c[NRCELLS-1].x-pow(sin(theta),2.)*p[i].r*p[i].r,0.5);
	 
	 	p[i].x= p[i].x+mfp*p[i].kx;
	 	p[i].y= p[i].y+mfp*p[i].ky;
   	 	p[i].z= p[i].z+mfp*p[i].kz;
     	p[i].r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);
	 	find_photon_angles(p,i,c);
	 
//		printf("r: %.f, R: %.f, cR: %.f\n",p[i].r/AU,R_MAX/AU,c[NRCELLS-1].x/AU );
	 	p[i].r=c[NRCELLS-1].x;
     
	 	p[i].flag=1;	
	 
       }
     
   }

 if (index==-100)
  {     

     if (p[i].r>c[NRCELLS-1].x) {flag1=1;flag2=1;index=-1;}       
     else  
       
       {
	 
	 s=binary_search_cell(p[i].r, 0, NRCELLS-1, c);	    
	 
//	 printf("s=%d\n",s);

	 for (k=0;k<NLCELLS-1;k++)
	   
	   if (p[i].theta>=c[(int)(k*NRCELLS*NMCELLS)].y && p[i].theta<=c[(int)((k+1)*NRCELLS*NMCELLS)].y)  
	     {
	       flag1=1;break;
	     }	    
	 
	 for(j=0;j<NMCELLS-1;j++)
	   
		 if (p[i].phi>= c[(int)(NRCELLS*j)].z && p[i].phi<=c[(int)(NRCELLS*(j+1))].z)
	     {	       
	       flag2=1;break;
	     }	    	 
	 
	 index=s+(k*NMCELLS+j)*NRCELLS;      
       }    
     
  //	 printf("index=%d r=%.1f\n",index,c[(int)index].x/AU);
//	 exit(0);

     if (flag1==0 || flag2==0)    

       printf("find.cell warning 1: r:%.6e  theta:%.6e phi:%.8e rmac:%.6e i:%d j:%d k:%d f1:%d f2:%d \n",
	      p[i].r, p[i].theta, p[i].phi,R_MAX*1.,s,k,j,flag1,flag2); 
     
     if (p[i].flag==0 && p[i].r>c[NRCELLS-1].x) 
	 {
       printf("find.cell warning 2: r:%.16e  theta:%.6e phi:%.6e rmax:%.16e  crmax:%.16e i:%d j:%d k:%d\n",
       p[i].r, p[i].theta, p[i].phi,R_MAX*1.,c[NRCELLS-1].x, s,k,j);
  printf("%.26e  %.26e %.26e \n",c[NRCELLS-1].x,R_MAX,R_MAX_CORE);
}		

     p[i].flag=1;
  }
  
 
 return (index);
}


/******************** modified binary search **************************/

int binary_search_cell (r, low, high, c)


double r;
int low;
int high;
struct cell_type c[];

{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (r>=c[middle].x && r<=c[middle+1].x) 

    return middle;

 else if (r<c[middle].x)

         high=middle-1;

 else low=middle+1;

 }

 printf("r:%.7e c:%.7e R_MAX:%.7e \n",r, c[NRCELLS-1].x,R_MAX);
printf("Problem in binary search\n");
exit(0);
return -1; /* search key not found */

}


/**********  origin ********************/

void find_photon_origin(p,i,star,s)

 struct photon_packet_type p[];
 int i;
 struct star_type star[];
 unsigned short s;
{ 
  if (p[i].lr<R_ORIGIN) p[i].origin=1;
  else p[i].origin=0;

  return;
}

/***** find theta-phi ************************************************************/
void find_photon_angles(p,i,c)

     struct photon_packet_type p[];
     int i;
     struct cell_type c[];
{
  
  double aux;

  if (p[i].r==0.) {p[i].theta=0.;  p[i].phi=0.;}     

  else  

    {
      p[i].theta=acos(p[i].z/p[i].r);
      
	if (p[i].x==0) {p[i].phi=0.;}

      else

	{p[i].phi=atan2(p[i].y,p[i].x);}

      if  (p[i].phi<0) {p[i].phi=p[i].phi+2.*M_PI;}
   
    }

 
  /* to avoid rounding errors */

  if (p[i].phi>c[(int)(NRCELLS*(NMCELLS-1))].z) p[i].phi=c[(int)(NRCELLS*(NMCELLS-1))].z  ;
  if (p[i].theta>c[(int)((NLCELLS-1)*NRCELLS*NMCELLS)].y) p[i].theta=c[(int)((NLCELLS-1)*NRCELLS*NMCELLS)].y;

  /* printf("%.10e %.10e\n", c[(int)(NRCELLS*(NMCELLS-1))].z,c[(int)((NLCELLS-1)*NRCELLS*NMCELLS)].y); */
  
  return;
}	  
