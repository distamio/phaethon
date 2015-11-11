#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "params.grid.h"
#include "params.h"
#include "constants.h"
#include "inout.h"
#include "opacity.h"
#include "phaethon.functions.h"
#include "grid.functions.h"


double MIN(double, double, double, double);


/***** function advance_photon **************************************/

/* advances photon till tau=0 or till it goes outside the grid */
/* and calculates the new cell (index) of the photon           */

void advance_photon(p,i,c,ncells,star,s)
     struct photon_packet_type p[];
     int i;
     struct cell_type c[];
     int ncells;
     struct star_type star[];
     unsigned short s;

{
  
  double mfp,density,b;
  int index;
  double total_mfp=0,lambda;
  double dtau=0,aux;
  double r_min,x_min,y_min,z_min,min_mfp,aux_max;
  double stepsize=0.,theta,ltau,theta_be,phi,omega,d,dens_step;

  if (c[p[i].cell].mass==0) 
   {
    if (p[i].r<R_DUST)  {advance_photon_free(p,i,c,ncells,0.5*(R_DUST-p[i].r)+0.5*(c[2].x-c[1].x),star,s);}
    else
    advance_photon_free(p,i,c,ncells,0.5*(c[2].x-c[1].x),star,s);
   }

  /* check minimum distance photon travels for max dens and for max absorption */

  b=pow(p[i].r*p[i].r- pow(p[i].x*p[i].kx+p[i].y*p[i].ky+p[i].z*p[i].kz,2.),0.5);

  if (b<R_DUST) b=c[1].x;

  
#ifndef DUST_MIX
  aux_max=(p[i].k_abs+p[i].s_scat)*dens(b,0.5*M_PI,0); // this will not work for comet like!!!!!
  if (DIM==4) {printf("\n:: CHECK advance.c :: FF9\n"); exit(0);}
#else
  if ((p[i].k_abs+p[i].s_scat)>(p[i].k_abs2+p[i].s_scat2)) 
	
 	aux_max=(p[i].k_abs+p[i].s_scat)*dens(b,0.5*M_PI,0);
  else
	
	
 	aux_max=(p[i].k_abs2+p[i].s_scat2)*dens(b,0.5*M_PI,0);
#endif

  min_mfp=p[i].tau/aux_max;

  x_min=p[i].x+p[i].kx*min_mfp;
  y_min=p[i].y+p[i].ky*min_mfp;
  z_min=p[i].z+p[i].kz*min_mfp;

  r_min=pow(x_min*x_min+y_min*y_min+z_min*z_min,0.5);


  if (r_min>c[NRCELLS-1].x) p[i].cell=-1;

 
 else
 
 {
     find_photon_angles(p,i,c);

#ifdef DUST_MIX
     find_photon_dust_kind(p,i);
#endif     

      while(p[i].cell>0 && p[i].tau>0)
	{
	  density=dens(p[i].r, p[i].theta, p[i].phi); 

#ifdef DUST_MIX

          if (p[i].dust_kind==1)
{
            aux=(p[i].k_abs2+p[i].s_scat2)*density;
}
          else
#endif
{
	  aux=(p[i].k_abs+p[i].s_scat)*density;    
}
	  
         stepsize=c[p[i].cell+1].x-c[p[i].cell].x;

         if (stepsize<0) stepsize=c[p[i].cell].x-c[p[i].cell-1].x;
 
	  mfp=MIN(STEP_CELL*stepsize,STEP_R*p[i].r,STEP_MFP/aux,(p[i].tau+1e-9)/aux);

	  dtau=mfp*aux;          

	  p[i].x=p[i].x+p[i].kx*mfp;
	  p[i].y=p[i].y+p[i].ky*mfp;
	  p[i].z=p[i].z+p[i].kz*mfp;
	  p[i].r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);

	  find_photon_angles(p,i,c);

#ifdef DUST_MIX
          find_photon_dust_kind(p,i);
#endif  
	  p[i].tau=p[i].tau-dtau;     
	
          p[i].cell=find_photon_cell(p,i,c,star,s);
		  		  

 if (c[p[i].cell].mass==0)
   {
    if (p[i].r<R_DUST)  advance_photon_free(p,i,c,ncells,0.5*(R_DUST-p[i].r)+0.5*(c[2].x-c[1].x),star,s);
    else
    advance_photon_free(p,i,c,ncells,0.5*(c[2].x-c[1].x),star,s);
   }

      }

}
  return;
}


/***** function advance_photon_free **************************************/

/* advances photon till it meets a cell with non-zero density and returns
   the cell index of the new cell */

/* no alteration in tau           */

/* returns -1 if the photon goes outside the grid */

void advance_photon_free(p,i,c,ncells,mfp, star,s)

 struct photon_packet_type p[];
 int i;
 struct cell_type c[];
 int ncells;
 double mfp;
 struct star_type star[];
 unsigned short    s;

{
  while (c[p[i].cell].mass==0 && (p[i].cell>=0)) 
        {     
// printf("%.4e %.4e here\n", p[i].r/AU,mfp/AU);
          p[i].x=p[i].x+p[i].kx*mfp;
          p[i].y=p[i].y+p[i].ky*mfp;
          p[i].z=p[i].z+p[i].kz*mfp;         
          p[i].r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);
          p[i].cell=find_photon_cell(p,i,c,star,s);          
         }    
 return;
}



/***** function find_mfp_free  **************************************/

double find_mfp_free(c,ncells)

 struct cell_type c[];

 int ncells;

{
  int i;
  double mfp=c[1].x-c[0].x; 

  for (i=1;i<ncells-1;i++)

       if (mfp>c[i+1].x-c[i].x) mfp=c[i+1].x-c[i].x;

  mfp=mfp*STEP_CELL;

 return (mfp);
}
/***** function find_step_size  **************************************/

double find_stepsize(p,i,c)

 struct photon_packet_type p[];
 long i;
 struct cell_type c[];

{
 double stepsize;

       /* calculate step size according to the width of the local cells */
  
     if (p[i].cell==0)   

	 {
          if ((c[1].x-c[0].x)<(c[2].x-c[1].x)) stepsize=STEP_CELL*(c[1].x-c[0].x);

          else stepsize=STEP_CELL*(c[2].x-c[1].x);
	 }

       else 
         { 
          if ((c[p[i].cell].x-c[p[i].cell-1].x)<(c[p[i].cell+1].x-c[p[i].cell].x)) 

	    
          stepsize=STEP_CELL*(c[p[i].cell].x-c[p[i].cell-1].x);

          else stepsize=STEP_CELL*(c[p[i].cell+1].x-c[p[i].cell].x);
	           
	 }

     return (stepsize);


}


/***** function advance_photon **************************************/

/* advances photon till tau=0 or till it goes outside the grid */
/* and calculates the new cell (index) of the photon           */

void advance_photon2(p,i,c,ncells,mfp_free,star, s)

 struct photon_packet_type p[];
 long i;
 struct cell_type c[];
 double mfp_free;
 int ncells;
 struct star_type star[];
  unsigned short s;

{
 double mfp,density,r;
 int index;
 double lambda;
 double dtau=0;
 double x0,y0,z0,x1,y1,z1,density1,alpha,r1;
 double r_min,x_min,y_min,z_min,min_mfp,aux_max,dtau_big_step,BIG_STEP;
 int cell1;
  double stepsize=0.;
 double moment[10];

     if (c[p[i].cell].mass==0)
       
       advance_photon_free(p,i,c,ncells,mfp_free,star,s);                      

/* calculate */

       aux_max=(p[i].k_abs+p[i].s_scat)*c[0].dens; 
       min_mfp=p[i].tau/aux_max; 

       x_min=p[i].x+p[i].kx*min_mfp; 
       y_min=p[i].y+p[i].ky*min_mfp; 
       z_min=p[i].z+p[i].kz*min_mfp; 
 
       r_min=pow(x_min*x_min+y_min*y_min+z_min*z_min,0.5); 
       
       if (r_min>R_MAX) p[i].cell=-1; 

       else

	 {
	   printf("ok\n");

          r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);

          /* linear interpolation for dens/ needs the last cell to have a density!!! */ 

       density=((r-c[p[i].cell].x)*c[p[i].cell+1].dens+  
             (c[p[i].cell+1].x-r)*c[p[i].cell].dens)/  
	     (c[p[i].cell+1].x-c[p[i].cell].x);  

      x0=p[i].x;
      y0=p[i].y;
      z0=p[i].z;

      dtau_big_step=0.;
      mfp=0;
       
   
      
     while (p[i].tau>dtau_big_step){

      BIG_STEP=5*find_stepsize(p,i,c);
      x1=x0+p[i].kx*BIG_STEP; 
      y1=y0+p[i].ky*BIG_STEP; 
      z1=z0+p[i].kz*BIG_STEP; 
      
      r1=pow(x1*x1+y1*y1+z1*z1,0.5); 

      cell1=find_cell_from_position(r1,c);
       
      if (cell1<0)
	{
       cell1=NCELLS-2;
       r1=R_MAX;
	}

      density1=((r1-c[cell1].x)*c[cell1+1].dens+  
             (c[cell1+1].x-r1)*c[cell1].dens)/  
	     (c[cell1+1].x-c[cell1].x);  

      dtau_big_step=(p[i].k_abs+p[i].s_scat)*(density+density1)*(BIG_STEP)/2.;
 
      if (p[i].tau>dtau_big_step)
	{
        p[i].tau=p[i].tau-dtau_big_step;       
        mfp=mfp+BIG_STEP;
        density=density1;
        x0=x1;
        y0=y1;
        z0=z1;
	}
     }

      alpha=(density1-density)/(BIG_STEP);
      mfp=mfp+(-density/alpha)+pow(density*density/(alpha*alpha)+2.*p[i].tau/((p[i].k_abs+p[i].s_scat)*alpha),0.5);
      dtau=p[i].tau;

       p[i].x=p[i].x+p[i].kx*mfp;
       p[i].y=p[i].y+p[i].ky*mfp;
       p[i].z=p[i].z+p[i].kz*mfp;

       
       p[i].tau=p[i].tau-dtau; 

       p[i].cell=find_photon_cell(p,i,c,star,s);

       /*
   if (p[i].cell<0) break; 

   if (p[i].tau<=0 && c[p[i].cell].dens>0) break;                   
 
   if (c[p[i].cell].dens==0) advance_photon_free(p,i,c,ncells,mfp_free);
       */
	 }
	

 return;
}

/**************************************************************************************************/
int find_cell_from_position(r,c)

     double r;
     struct cell_type c[];

{
 int s,flag=0;
 int index=-100;

/* 1D spherical grid */

if (r>R_MAX) return (-1);

else

for (s=NCELLS-1;s>0 ;s--)
      
if ( r>c[s-1].x && r<=c[s].x ) {index=s-1;break;}                	       
     
return (index);
   
}


/********************* minimum ************************************************/

double MIN(x,y,z,w)

     double x,y,z,w;

{
  double min;

  if (x>y) min=y;
  else min=x;

  if (z<min) min=z;

  if (w<min) min=w;

  return (min);

}



 	    

/****** ADVANCE SHADOW PHOTON ************************************/	    

void advance_shadow_photon(p,i,c,ncells,mfp_free)
     
     struct photon_packet_type p[];
     long i;
     struct cell_type c[];
     double ncells;
     double mfp_free;
     
{
  
  double mfp,density;
  int index;
  double total_mfp=0,lambda;
  double dtau=0,aux;
  double r_min,x_min,y_min,z_min,min_mfp,aux_max;
  double stepsize=0.,theta,ltau,theta_be,phi,omega,d,dens_step;
  
 

#ifdef INC_SHADOW_PHOTONS 
    p[i].tau=0.;

    /* printf("%.4e %d %.4e\n",p[i].r,p[i].cell,p[i].tau); */
    while(p[i].r<=R_MAX && (p[i].cell)>=0 && p[i].tau<SHADOW_TAU_MAX)

      {
	/*	printf("%.4e %.4e %.4e\n",p[i].r,p[i].theta,p[i].phi); */
	
	density=dens(p[i].r, p[i].theta, p[i].phi);
	aux=(p[i].k_abs+p[i].s_scat)*density;
	
	dens_step=p[i].r/2;
	
	/*	printf("%.4e %.4e %.4e %.4e\n",STEP_DENS*dens_step,STEP_R*p[i].r,STEP_MFP/aux,(p[i].tau+1e-9)/aux); */	
	
	mfp=MIN(STEP_DENS*dens_step,STEP_R*p[i].r,STEP_MFP/aux,STEP_MFP/aux);
	/* printf("-----%.4e\n",mfp); */
    
	if (mfp==0) mfp=0.1*(c[p[i].cell+1].x-c[p[i].cell].x); /*** BE CAREFUL **/

	dtau=mfp*aux;          
	/* printf("r_bef %.4e\n",p[i].r); */
	p[i].x=p[i].x+p[i].kx*mfp;
	p[i].y=p[i].y+p[i].ky*mfp;
	p[i].z=p[i].z+p[i].kz*mfp;
	p[i].r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);
	/* printf("kx ky kz in %.4e %.4e %.4e\n",p[i].kx,p[i].ky,p[i].kz); */
	p[i].cell=find_photon_cell(p,i,c); 
	
	p[i].tau=p[i].tau+dtau;     
      }

#endif
	
      return;
}

