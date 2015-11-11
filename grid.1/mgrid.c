
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "structs.h"
#include "params.h"
#include "phaethon.params.h"
#include "grid.functions.h"
#include "inout.h"

void mgrid(c)

  struct cell_type  c[];
{
struct cell_type  ctype;
int s,i,k,j;
double mc_density; /* molecular cloud density */
double tau[NRCELLS+1];
double total_tau=0;

FILE *fp;


printf("---------- MAKING RT GRID ---------- X ---------------------\n");
if (DISC==1) printf (":: GRID WITH DISC\n");
if (ENVELOPE==1) printf (":: GRID WITH ENVELOPE\n");
if (ENVELOPE==0) printf (":: GRID WITHOUT ENVELOPE\n");
/* initialize cells' parameters */

if (ENVELOPE==0) {R_MAX_CORE=DISC_R;R_MAX=DISC_R;}

init_cell(c,R_DUST,R_MAX_CORE);

/* compute density  */

find_cell_dens(c);

/* compute cell mass */

find_cell_mass(c);

/* compute ambient molecular Cloud quantities */

 if (MCCELLS>0 && ENVELOPE==1) MC_fix_cells(c,R_MAX_CORE,R_MAX);

/* output the incells.dat file*/

write_cells(c,NCELLS);


/* routine to find A_V */

  
 if (( fp = fopen("Av.r.dat", "w")) == NULL)
   {
     fprintf(stderr, "FAILED TO OPEN Av.r.dat FILE. \n");
     exit(0);
   }
 
 tau[0]=0;

 
 for (i = 1; i<=NRCELLS; i++) 
   
   tau[i]=c[NRCELLS-i-1].tau+tau[i-1];
   
 fprintf(fp,"  A_v     r (AU)\n"); 
  
 for (i =0; i <=NRCELLS; i++) 
   { 
//   fprintf(fp,"A_v=%.3f  r=%.4e \n", tau[i]*1.086, c[NRCELLS-i-1].x/AU);
     fprintf(fp,"%.3f  %.4e \n", tau[i]*1.086, c[NRCELLS-i-1].x/AU);
  }

fclose(fp);
/* Append R_MAX at the end of grid.params.h file */

if (( fp = fopen("params.grid.h", "a+")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN params.grid.h FILE. \n");
         exit(0);
         }

fclose(fp);

 find_tau();

if (DISC==1) find_tau_disc();
}

/*************** initialize cell values ******************************/


void init_cell(c,rmin,rmax)

     struct cell_type c[];
     double rmin,rmax;
 
{
  int i,j,k;
  int s,m,l=0;
  double dr;
  double xtable[BECELLS];
  double ytable[NLCELLS];
  double ztable[NMCELLS];
  double dx,dy,dz;

  xtable[0]=0;
  xtable[1]=rmin;

if (LINEAR==1) {
  dx=(rmax-rmin)/(BECELLS-2.);
  for (i=2;i<BECELLS;i++)  xtable[i]=(i-1)*dx+rmin;     
  }
else
 {
   dx=log10(rmax/rmin)/(BECELLS-2.);
   for (i=2;i<BECELLS;i++) xtable[i]=pow(10,dx)*xtable[i-1];
 }
   xtable[BECELLS-1]=R_MAX_CORE; // to avoid stupid rounding errors 

  ytable[0]=0;

  /* iso-angle cells */

  if (ISOANGLES==1)
    {
      dy=M_PI/(NLCELLS-1);
      
      for (i=0;i<NLCELLS;i++) 
	{      	
	  ytable[i]=i*dy; 
	}
    }
  else 
    {
      /* iso-cos  cells */
      
      dy=2./(NLCELLS-1.);
      
      for (i=0;i<NLCELLS;i++) 
	{      	
	  ytable[i]=acos(1-i*dy); 
	} 
    }

  /* iso-phi cells */
  ztable[0]=0;
 
 dz=2*M_PI/(NMCELLS-1);
  for (i=0;i<NMCELLS;i++) 
    {      	
      ztable[i]=i*dz; 
    }




  for(s=0;s<NRCELLS;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
	{
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x=xtable[s];
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].y=ytable[k];
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].z=ztable[j];
	    
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].nabs=0.; 
	  
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].abslum=0; 
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].next_n=0; 
	}

return;
}

/************** find tau of the cell *****************************/

void find_cell_tau(c)

 struct cell_type c[];

{
 int s,i,k,j;
 double r;


 for(s=0;s<BECELLS+NCELLS-1;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
        {
          r=(c[s].x+c[s+1].x)/2;
          c[(int)(s+(k*NMCELLS+j)*NRCELLS)].tau=c[(int)(s+(k*NMCELLS+j)*NRCELLS)].dens*c[(int)(s+(k*NMCELLS+j)*NRCELLS)].K_VISUAL*(-c[s].x+c[s+1].x);
        }
return;
}
    

/************** find density of the cell *****************************/

void find_cell_dens(c)

 struct cell_type c[];

{
 int s,i,k,j;
 double r;
 double total_tau=0,r1,r2;

 for(s=0;s<BECELLS;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
	{
          r=c[s].x;
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].dens=
	    dens(c[s].x,c[(int)(k*NMCELLS*NRCELLS)].y,c[(int)(j*NRCELLS)].z);
	}

return;
}



/*********** find mass of each cell *****************************************/

void find_cell_mass(c)

 struct cell_type c[];

{
 int s,i,k,j;
 double total_mass=0.;
 double disc_mass=0.;



 for(s=0;s<BECELLS-1;s++)
    for(k=0;k<NLCELLS-1;k++)
      for(j=0;j<NMCELLS-1;j++)
	{
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass=
	    mass( c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x, c[(int)(s+1+(k*NMCELLS+j)*NRCELLS)].x,
		  (c[(int)(s+(k*NMCELLS+j)*NRCELLS)].y),     (c[(int)(s+((k+1)*NMCELLS+j)*NRCELLS)].y),
		  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].z, c[(int)(s+(k*NMCELLS+j+1)*NRCELLS)].z);
 
	  total_mass=total_mass+c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass;

if (DISC==1)  
{
if ( c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x<DISC_R)
           disc_mass=disc_mass+c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass;
}

	}
 



 printf("\nTotal mass (from cells) 		: %.3e\n",total_mass/SOLAR_MASS);

if (DISC==1) 
 printf(               "Total mass (from density profile)	: %.5e %.5e %.5e\n",
mass(0,R_MAX_CORE,0,0.5*M_PI-DISC_THETA_OPENING*DISC_ALPHA-1e-10,0,2*M_PI)/SOLAR_MASS,
mass(0,R_MAX_CORE,0.5*M_PI-DISC_THETA_OPENING*DISC_ALPHA,0.5*M_PI+DISC_THETA_OPENING*DISC_ALPHA,0,2*M_PI)/SOLAR_MASS,
mass(0,R_MAX_CORE,0.5*M_PI+DISC_THETA_OPENING*DISC_ALPHA+1e-10,M_PI,0,2*M_PI)/SOLAR_MASS);
else 
{
if (JET==1) 
 printf(               "Total mass (from density profile)	:  need to do it.... special care needed\n");
else

 printf(               "Total mass (from density profile)	: %.5e\n",
mass(0,R_MAX_CORE,0,M_PI,0,2*M_PI)/SOLAR_MASS);
}


if (DISC==1)  printf("\nTotal disc Mass (from cells)            : %.3e\n",disc_mass/SOLAR_MASS);



return;

}
			   

  
/*************** initialize cell values ******************************/


void MC_fix_cells(c,rmin,rmax) 

 struct cell_type c[]; 
 double rmin;
 double rmax;
 
{
  int i,j,k;
  int s,m,l=0;
  double dr;
  double xtable[NRCELLS];
  double ytable[NLCELLS];
  double ztable[NMCELLS];
  double dx,dy,dz,total_mass=0;;
  double total_tau=0;

  xtable[0]=0;

  dx=(rmax-rmin)/MCCELLS;

//  printf("\nrmin %.4e rmax %.4e dx %.4e\n",rmin,rmax,dx); 

  for (i=BECELLS-1;i<NRCELLS;i++) 
    {      	
      xtable[i]=(i-BECELLS+1.)*dx+rmin; 
 //    printf("\nrmin %.4e\n",xtable[i]); 
    }

  xtable[NRCELLS-1]=R_MAX; // to avoid stupid rounding errors

  ytable[0]=0;



  /* iso-angle cells */

  if (ISOANGLES==1)
    {
      dy=M_PI/(NLCELLS-1);

      for (i=0;i<NLCELLS;i++) 
	{      	
	  ytable[i]=i*dy; 
	}
    }

  else

  /* iso-cos  cells */
    {
      dy=2./(NLCELLS-1.);
      
      for (i=0;i<NLCELLS;i++) 
	{      	
	  ytable[i]=acos(1-i*dy); 
	} 
    }

  

  /* iso-phi cells */
  ztable[0]=0;
 
 dz=2*M_PI/(NMCELLS-1);
  for (i=0;i<NMCELLS;i++) 
    {      	
      ztable[i]=i*dz; 
    }


  for(s=BECELLS-1;s<NRCELLS;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
	{
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x=xtable[s];
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].y=ytable[k];
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].z=ztable[j];
	}
  
  for(s=BECELLS-1;s<NRCELLS;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
	{        
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].temp=0.; 
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].dens=
	    dens(c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x,c[(int)(s+(k*NMCELLS+j)*NRCELLS)].y,
		       c[(int)(s+(k*NMCELLS+j)*NRCELLS)].z); 
	}
  
  
  for (s=BECELLS-1;s<NRCELLS-1;s++) 
    for(k=0;k<NLCELLS-1;k++)
      for(j=0;j<NMCELLS-1;j++)
	{      	     
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass=0.;

	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass=
	    mass( c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x, c[(int)(s+1+(k*NMCELLS+j)*NRCELLS)].x,
		  (c[(int)(s+(k*NMCELLS+j)*NRCELLS)].y),     (c[(int)(s+((k+1)*NMCELLS+j)*NRCELLS)].y),
		  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].z, c[(int)(s+(k*NMCELLS+j+1)*NRCELLS)].z);
	  
	  total_mass=total_mass+c[(int)(s+(k*NMCELLS+j)*NRCELLS)].mass;
	  
//	  printf("%.4e %.4e\n", c[(int)(s+(k*NMCELLS+j)*NRCELLS)].x,c[(int)(s+1+(k*NMCELLS+j)*NRCELLS)].x);
	  
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].nabs=0.; 
	  
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].abslum=0; 
	  c[(int)(s+(k*NMCELLS+j)*NRCELLS)].next_n=0; 
	} 
  

 for (s=BECELLS-1;s<NRCELLS-1;s++)
    for(k=0;k<NLCELLS-1;k++)
      for(j=0;j<NMCELLS-1;j++)
        {
        c[(int)(s+(k*NMCELLS+j)*NRCELLS)].tau= 
        c[(int)(s+(k*NMCELLS+j)*NRCELLS)].K_VISUAL*
        c[(int)(s+(k*NMCELLS+j)*NRCELLS)].dens*
       (c[s+1].x-c[s].x); 
        }


for (s=BECELLS-1;s<NRCELLS-1;s++)
        total_tau=total_tau+c[s].tau;



 printf("\nTotal cloud vis tau: %.3f\n",total_tau);
  printf("\nTotal cloud Mass  : %.3e\n",total_mass/SOLAR_MASS);
  printf("Total cloud Mass 2: %.5e\n",mass(R_MAX_CORE,R_MAX,0,M_PI,0,2*M_PI)/SOLAR_MASS); 
return;
}







/******* FIND TAU ********************************************************/
void find_tau()

    
{
 double dtau=0.1;
 double r=0., tau=0.,tau2=0.,sum=0. ,mfp=0.,tau_eq,tau_eq2,tau_eq2_ext;
 double tau_pole,tau_pole2,tau_pole2_ext;
 double tau_spole,tau_spole2,tau_spole2_ext;
 double tau_R0=0.,tau_ext=0.;
 double theta=0.,k_vis;
 printf("\n:: tau visual (at different directions) \n");
 for (theta=0;theta<=180;theta=theta+10)
   {
    
     r=R_DUST;


     while (r<R_MAX_CORE)

       {	 
#ifndef DUST_MIX
     k_vis=K_VISUAL0;
#else
     k_vis=find_k_visual(r);
//   printf("k_vis=%.4e\n",k_vis);

#endif
	 mfp=dtau/(k_vis*dens(r,(M_PI*theta/180.),0));	 
	 tau=tau+dtau;
         if (r>R0) tau_ext=tau_ext+dtau;
         if (r<R0) tau_R0=tau_R0+dtau;
	 r=r+mfp;	 
	// printf("tau=%.4f r=%.3f \n",find_k_visual(r),r/AU);
	 if (theta==0) tau_pole2=tau; 
	 if (theta==90) tau_eq2=tau;   
	 if (theta==180) tau_spole2=tau; 
	 if (theta==0) tau_pole2_ext=tau_ext; 
	 if (theta==90) tau_eq2_ext=tau_ext;   
	 if (theta==180) tau_spole2_ext=tau_ext; 
       }     


     while (r<R_MAX)

       {
#ifndef DUST_MIX
     k_vis=K_VISUAL0;
#else
     k_vis=find_k_visual(r);
#endif

       	 mfp=dtau/(k_vis*dens(r,(M_PI*theta/180.),0));       
	 tau2=tau2+dtau;
	 r=r+mfp;
	 /* printf("tau=%.4f r=%.3f density= %.3e\n",tau,r/R_DUST,dens(r,theta)); */
       }
    
    printf(" theta=%3.f :: tau_R0=%8.4f  tau_core=%8.4f   tau_cloud=%8.4f     \n", theta, tau_R0, tau, tau2);

     
     if (theta==0) tau_pole=tau2; 
     if (theta==90) tau_eq=tau2; 
     if (theta==180) tau_spole=tau2; 
     tau=0;tau2=0;tau_ext=0;tau_R0=0;

   }

 printf("\n"); 
 printf(" e_cloud(eq/pole)    =%6.3f   e(eq/pole)   = %6.3f  \n",tau_eq/tau_pole,(tau_eq+tau_eq2)/(tau_pole+tau_pole2));  
 if (theta>90) 
   printf(" e_cloud(spole/pole) =%6.3f   e(spole/pole)= %6.3f  \n",
	  tau_spole/tau_pole,(tau_spole+tau_spole2)/(tau_pole+tau_pole2)); 
 printf("\n"); 

 
printf(" *Exluding* region within R0 and outside R_MAX_CORE\n"); 
    printf("  e(eq/pole)   = %6.3f \n",  tau_eq2_ext/tau_pole2_ext);  
 if (theta>90) 
    printf("  e(spole/pole)= %6.3f  \n", tau_spole2_ext/tau_pole2_ext); 
 printf("\n"); 

 return;
 
 
}


/******* FIND TAU DISC ********************************************************/
void find_tau_disc()


{
 double dtau=0.1;
 double r=R_DUST, tau=0.,tau2=0.,sum=0. ,mfp=0.,tau_eq,tau_eq2;
 double tau_pole,tau_pole2;
 double tau_spole,tau_spole2;
 double theta=0.,k_vis;
 printf("\n-----------------tau visual----------------------\n");
 for (theta=0;theta<=180;theta=theta+30)
   {
     r=R_DUST;

     while (r<DISC_R)

       {

#ifndef DUST_MIX
     k_vis=K_VISUAL0;
#else
     k_vis=find_k_visual(r);
#endif
         mfp=dtau/(k_vis*dens(r,M_PI*theta/180.,0));
         tau=tau+dtau;
         r=r+mfp;
         /* printf("tau=%.4f r=%.3f density= %.3e\n",tau,r/R_DUST,dens(r,theta,0));*/
         if (theta==0) tau_pole2=tau;
         if (theta==90) tau_eq2=tau;
         if (theta==180) tau_spole2=tau;
       }

     while (r<R_MAX)

       {
#ifndef DUST_MIX
     k_vis=K_VISUAL0;
#else
     k_vis=find_k_visual(r);
#endif
         mfp=dtau/(k_vis*dens(r,(M_PI*theta/180.),0));
         tau2=tau2+dtau;
         r=r+mfp;
         /* printf("tau=%.4f r=%.3f density= %.3e\n",tau,r/R_DUST,dens(r,theta)); */
       }

     printf(" (theta=%3.f) tau_disc=%10.4f     \n",theta, tau, tau2);


     if (theta==0) tau_pole=tau2;
     if (theta==90) tau_eq=tau2; 
     if (theta==180) tau_spole=tau2;

     tau=0;tau2=0.;
   }

  
 return;

 
}

