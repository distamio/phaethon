
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#include "structs.h"
#include "params.grid.h"
#include "params.opa.h"
#include "constants.h"
#include "params.h"
#include "phaethon.params.h"

#include "functions.h"
#include "inout.h"
#include "opacity.h"
#include "phaethon.functions.h"
#include "photons.h"
#include "spec.functions.h"
#include "interactions.h"

// *** find the number of photons needed to be absorbed for the next temperature increase  ***

#ifndef DUST_MIX
void find_cells_next_n(c, ncells, lum, mean_op_ptr)


 struct cell_type c[];
 double lum;
 long ncells;
 struct mean_op_type *mean_op_ptr;

 {
  int cell;

    for (cell=0;cell<ncells;cell++)

      c[cell].next_n=invree(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op_ptr);

return;
}
#else

void find_cells_next_n(c, ncells, lum, mean_op_ptr,mean_op2_ptr)


 struct cell_type c[];
 long ncells;
 double lum;
 struct mean_op_type *mean_op_ptr;
 struct mean_op_type *mean_op2_ptr;
 {
  long cell;

    for (cell=0;cell<ncells;cell++)
	{	
if (cell_dustkind(c,cell)==0)

      c[cell].next_n=invree(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op_ptr);

else
      c[cell].next_n=invree2(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op2_ptr);
}

return;
}

#endif

/************************ function absorb_photon **************************************/
double absorb_photon(lum,cell, c,mean_op_ptr)


 struct cell_type c[];
 struct mean_op_type * mean_op_ptr;
 long cell;
 double lum; 

{
 double old_temp; /* temperature before absorbing the photon */

#ifdef INFO_OUTPUT   
 printf(" ABSORBING...\n");
#endif

 if (c[cell].mass==0) 
   { 
      printf("Mass of absorbing cell is zero !!");
      printf("%ld\n",cell);
      exit(0);
    }
 old_temp=c[cell].temp; 

 c[cell].abslum=c[cell].abslum+lum;

 /* count how may photons have been absorbed */
 c[cell].nabs++;

/* find new temperature after absorbing the photon packet */
#ifdef INFO_OUTPUT
 printf("old Temp: %.1f [cell: %d]\n",c[cell].temp,cell);
#endif

 if (c[cell].nabs==1)
   {
     c[cell].temp=solve_ree(c[cell].abslum,c[cell].mass,old_temp,c[cell].nabs-1,mean_op_ptr);
     c[cell].next_n=invree(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op_ptr);

     /*printf("%d %f \n", c[cell].next_n,(MAX_TEMP-1.)/(TEMPDIM-1.) ); */
   }
 else
   {
     if (c[cell].nabs==c[cell].next_n && c[cell].temp<MAX_TEMP-1)
       { /* printf("+++++++++++++++++++++++++%d %d\n", c[p[i].cell].nabs,c[p[i].cell].next_n );  */ 
	 c[cell].temp=c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.);
	 c[cell].next_n=invree(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op_ptr);

	 /* printf("-------------------------%d %d\n", c[p[i].cell].nabs,c[p[i].cell].next_n ); */
       }
      else

       if (c[cell].nabs>c[cell].next_n && c[cell].temp<MAX_TEMP-1) 
           /* for the camse that one absorption leads to temp>t+DTmin */
	 {
	   c[cell].temp=solve_ree(c[cell].abslum,c[cell].mass,old_temp,c[cell].nabs-1,mean_op_ptr);
	   c[cell].next_n=invree(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op_ptr);
	 }
   }


#ifdef INFO_OUTPUT
 printf("new Temp: %.1f [cell: %d]\n",c[cell].temp,cell);
#endif

 return(c[cell].temp) ;
}

/************************ function photon_reemit  **************************************/

void reemit_photon(p,i,reemit_ptr)

 struct photon_packet_type p[];
 long i;
 struct reemit_type *reemit_ptr;
 
{
 double R2,R3,R4,aux;
 int index_aux=0;

#ifdef INFO_OUTPUT
 printf(" REEMITTING...\n");
 printf("old lambda (micron): %7.1f fpoint:%5d \n",1e4*c_CONST/p[i].f,p[i].fpoint);
#endif

/* calculate new frequency of the reemitted photon */
 find_new_f(p,i,reemit_ptr);   

 /*
    DOES NOT WORK!

    index_aux=(int)(p[i].temp*(MAX_TEMP-1)/(TEMPDIM-1));
    p[i].f=(*reemit_ptr).rfreq[index_aux]; 
    
    p[i].fpoint=(*reemit_ptr).fpoint[index_aux]; */

#ifdef INFO_OUTPUT
 printf("new lambda (micron): %7.1f fpoint:%5d\n",1e4*c_CONST/p[i].f,p[i].fpoint);
#endif
 
 /* give photon packet new direction and new tau */

 R2=(double) (rand())/RAND_MAX;
 R3=(double) (rand())/RAND_MAX;  
 R4=(double) (rand())/RAND_MAX;  

 /* printf( "%.4e %.4e %.4e %.4e\n",R2,R3,R4,1.*RAND_MAX );*/
 
 p[i].kz=2.*R2-1;

 aux=sin(acos(p[i].kz));

 p[i].kx=aux*cos(2.*M_PI*R3);
 p[i].ky=aux*sin(2.*M_PI*R3);

 p[i].tau=-log(R4);
 p[i].abs=p[i].abs+1; 

 /* p[i].abs=1; p[i].scat=0; */

#ifdef INFO_OUTPUT
 printf("tau: %.4f \n",p[i].tau);
#endif

 return;
}

/******** Function ree (radiative equllibrium equation) *************/


double ree(temp,abslum,mass,mean_op_ptr)

     double temp; 
     double abslum;
     double mass;
     struct mean_op_type * mean_op_ptr;
    

{
return 
(sigma_CONST*pow(temp,4)-(abslum/(4.*mean_opacity(temp,mean_op_ptr)*mass)));
}

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/***************************** Function invree ***************************/
/* returns the number of photons needed to be absorbed for the tempetarute
   of a cell with mass 'mass' to be raised to 'temp'
*/
 int invree(temp,dL,mass,mean_op_ptr)

     double temp; 
     double dL;
     double mass;
     struct mean_op_type * mean_op_ptr;

{
return ((int)(4*sigma_CONST*pow(temp,4.)*mean_opacity(temp,mean_op_ptr)*mass/dL));
}
  



/******* Function solve ree (Radiative Equilibrium Equation) ********/

double solve_ree(abslum,mass,old_temp,nabs,mean_op_ptr)

double abslum;
double mass;
double old_temp;
int nabs; /* previously absorbed photons, ie c.nabs-1 */ 
struct mean_op_type * mean_op_ptr;

{
double temp,aux;
 
 aux=(MAX_TEMP-1.)/(TEMPDIM-1.);
 
 
 if (old_temp==0) old_temp=2*aux;
 
 for (temp=old_temp;temp<MAX_TEMP-1;temp=temp+aux)
   
   if (ree(temp,abslum,mass,mean_op_ptr)>=0) break;
 
 return(temp);
 
}


/***************** function find_new_f *************************************/

/* finds the frequency of the reemitted photon packet */

void find_new_f(p,i,reemit_ptr)
  
 struct photon_packet_type p[];
 long   i;
 struct reemit_type * reemit_ptr;

{
 double random;
 int index;
 int low=0;
 int  high;  
 int middle;

 /* random=generate_random(&idum); */

random=(double)(rand())/RAND_MAX;

index=(int)(p[i].temp/((MAX_TEMP-1)/(TEMPDIM-1)));
high=(*reemit_ptr).last[index];
/* 
for(k=NFREQ-1;k>=0;k=k-1)

if (reemit.prob[index][k]<random) break;

p[i].f=reemit.freq[k];

return;
*/

 while (low<=high) {

 middle= (low+high)/2;

 if (random>=(*reemit_ptr).prob[index][middle] && random<(*reemit_ptr).prob[index][middle+1]) 

    break;

 else if (random<(*reemit_ptr).prob[index][middle])

         high=middle-1;

 else low=middle+1;

 }

p[i].f=(*reemit_ptr).freq[middle];
p[i].fpoint=middle;

return;

}



#ifdef DUST_MIX

/************************ function absorb_photon_reemit  **************************************/
double absorb_photon2(lum,cell, c,mean_op2_ptr)


 struct cell_type c[];
 struct mean_op_type * mean_op2_ptr;
 long cell;
 double lum;

{
 double old_temp; /* temperature before absorbing the photon */

#ifdef INFO_OUTPUT   
 printf(" ABSORBING 2...\n");
#endif

 if (c[cell].mass==0)
   {
      printf("Mass of absorbing cell is zero !!");
      printf("%ld\n",cell);
      exit(0);
    }
 old_temp=c[cell].temp;

 c[cell].abslum=c[cell].abslum+lum;

 /* count how may photons have been absorbed */
 c[cell].nabs++;

/* find new temperature after absorbing the photon packet */
#ifdef INFO_OUTPUT
 printf("old Temp: %.1f [cell: %d]\n",c[cell].temp,cell);
#endif

 if (c[cell].nabs==1)

   {
     c[cell].temp=solve_ree2(c[cell].abslum,c[cell].mass,old_temp,c[cell].nabs-1,mean_op2_ptr);
     c[cell].next_n=invree2(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op2_ptr);
     /*printf("%d %f \n", c[p[i].cell].next_n,(MAX_TEMP-1.)/(TEMPDIM-1.) );*/
   }
 else
   {
     if (c[cell].nabs==c[cell].next_n && c[cell].temp<MAX_TEMP-1)
       { /* printf("+++++++++++++++++++++++++%d %d\n", c[p[i].cell].nabs,c[p[i].cell].next_n );  */
         c[cell].temp=c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.);
         c[cell].next_n=invree2(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op2_ptr);
         /* printf("-------------------------%d %d\n", c[p[i].cell].nabs,c[p[i].cell].next_n ); */
       }
      else

       if (c[cell].nabs>c[cell].next_n && c[cell].temp<MAX_TEMP-1)
           /* for the case that one absorption leads to temp>t+DTmin */
         {
           c[cell].temp=solve_ree2(c[cell].abslum,c[cell].mass,old_temp,c[cell].nabs-1,mean_op2_ptr);
           c[cell].next_n=invree2(c[cell].temp+(MAX_TEMP-1.)/(TEMPDIM-1.),lum,c[cell].mass,mean_op2_ptr);
         }
   }


#ifdef INFO_OUTPUT
 printf("new Temp: %.1f [cell: %d]\n",c[cell].temp,cell);
#endif

 return(c[cell].temp) ;
}

/************************ function photon_reemit  **************************************/

void reemit_photon2(p,i,reemit2_ptr)

 struct photon_packet_type p[];
 long i;
 struct reemit_type *reemit2_ptr;

{
 double R2,R3,R4,aux;
 int index_aux=0;

#ifdef INFO_OUTPUT
 printf(" REEMITTING 2...\n");
 printf("old lambda (micron): %7.1f fpoint:%5d \n",1e4*c_CONST/p[i].f,p[i].fpoint);
#endif

/* calculate new frequency of the reemitted photon */
 find_new_f2(p,i,reemit2_ptr);

 /*

    DOES NOT WORK!

    index_aux=(int)(p[i].temp*(MAX_TEMP-1)/(TEMPDIM-1));
    p[i].f=(*reemit_ptr).rfreq[index_aux]; 
    
    p[i].fpoint=(*reemit_ptr).fpoint[index_aux]; */

#ifdef INFO_OUTPUT
 printf("new lambda (micron): %7.1f fpoint:%5d\n",1e4*c_CONST/p[i].f,p[i].fpoint);
#endif

 /* give photon packet new direction and new tau */

 R2=(double) (rand())/RAND_MAX;
 R3=(double) (rand())/RAND_MAX;
 R4=(double) (rand())/RAND_MAX;

 /* printf( "%.4e %.4e %.4e %.4e\n",R2,R3,R4,1.*RAND_MAX );*/

 p[i].kz=2.*R2-1;

 aux=sin(acos(p[i].kz));

 p[i].kx=aux*cos(2.*M_PI*R3);
 p[i].ky=aux*sin(2.*M_PI*R3);

 p[i].tau=-log(R4);
 p[i].abs=p[i].abs+1;

 /* p[i].abs=1; p[i].scat=0; */

#ifdef INFO_OUTPUT
 printf("tau: %.4f \n",p[i].tau);
#endif

 return;
}



/******** Function ree2 (radiative equllibrium equation) *************/


double ree2(temp,abslum,mass,mean_op2_ptr)

     double temp; 
     double abslum;
     double mass;
     struct mean_op_type *mean_op2_ptr;

{
return 
(sigma_CONST*pow(temp,4)-(abslum/(4.*mean_opacity(temp,mean_op2_ptr)*mass)));
}

 
  

/***************************** Function invree2 ***************************/
/* returns the number of photons needed to be absorbed for the tempetarute
   of a cell with mass 'mass' to be raised to 'temp'
*/
 int invree2(temp,dL,mass,mean_op2_ptr)

     double temp; 
     double dL;
     double mass;
    struct mean_op_type * mean_op2_ptr;

{
return ((int)(4*sigma_CONST*pow(temp,4.)*mean_opacity(temp,mean_op2_ptr)*mass/dL));
}
  


/******* Function solve ree (Radiative Equilibrium Equation) ********/

double solve_ree2(abslum,mass,old_temp,nabs,mean_op2_ptr)

double abslum;
double mass;
double old_temp;
int nabs; /* previously absorbed photons, ie c.nabs-1 */ 
struct mean_op_type *mean_op2_ptr;
{

double temp,result,aux;
   
 
 aux=(MAX_TEMP-1.)/(TEMPDIM-1.);
 
 
 if (old_temp==0) old_temp=2*aux;
 
 for (temp=old_temp;temp<MAX_TEMP-1;temp=temp+aux)
   
   if (ree2(temp,abslum,mass,mean_op2_ptr)>=0) break;
 
 return(temp);
 


}


/***************** function find_new_f *************************************/

/* finds the frequency of the reemitted photon packet */

void find_new_f2(p,i,reemit2_ptr)

 struct photon_packet_type p[];
 long   i;
 struct reemit_type * reemit2_ptr;

{
 double random;
 int index;
 int low=0;
 int  high;
 int middle;

 /* random=generate_random(&idum); */

random=(double)(rand())/RAND_MAX;

index=(int)(p[i].temp/((MAX_TEMP-1)/(TEMPDIM-1)));
high=(*reemit2_ptr).last[index];
/* 
for(k=NFREQ-1;k>=0;k=k-1)

if (reemit.prob[index][k]<random) break;

p[i].f=reemit.freq[k];

return;
*/



 while (low<=high) {

 middle= (low+high)/2;

 if (random>=(*reemit2_ptr).prob[index][middle] && random<(*reemit2_ptr).prob[index][middle+1]) 

    break;

 else if (random<(*reemit2_ptr).prob[index][middle])

         high=middle-1;

 else low=middle+1;

 }

p[i].f=(*reemit2_ptr).freq[middle];
p[i].fpoint=middle;

return;

}

#endif

/***** function scatter_photon ********************************************/

void scatter_photon(p,i)

 struct photon_packet_type p[];
 long i;

{
 double R2,R3,R4;
 double costheta=0;     /* scattering angle cosine */
 double phi;
 double ks_photon,kx_photon,ky_photon,ks,ky,sin_phi,cos_phi;
 double g2;

#ifdef INFO_OUTPUT
 printf(" SCATTERING...\n");
#endif

 R2=(double) (rand())/RAND_MAX;
 R3=(double) (rand())/RAND_MAX;
 R4=(double) (rand())/RAND_MAX;


 g2=G_SCATTER*G_SCATTER;

 if (G_SCATTER!=0.) 
 {
 costheta=(1.+g2-pow((1.-g2)/(1+G_SCATTER-2.*G_SCATTER*R2),2))/(2.*G_SCATTER);

 if (costheta<1.e-10) costheta=0.; /* scattering theta */
 }

 else

 costheta=2.*R2-1;

 /* costheta= new kz at photon's system */

 phi=2.*M_PI*R3;                          /* scattering phi */

 ks_photon=sin(acos(costheta));/* sin(theta_scattering) -direction*/
 
 /* new direction in photon's system */

 kx_photon=ks_photon*cos(phi);
 ky_photon=ks_photon*sin(phi);

 /* sin(theta_before_scattering) */

 ks=sin(acos(p[i].kz));

 /* new direction in original system */

 
 ky=p[i].kz*ky_photon+ks*costheta;
 p[i].kz=-ks*ky_photon+p[i].kz*costheta;

 sin_phi=(p[i].ky/ks);
 cos_phi=(p[i].kx/ks);

 p[i].kx=sin_phi*kx_photon+cos_phi*ky;
 p[i].ky=-cos_phi*kx_photon+sin_phi*ky;

 p[i].tau=-log(R4);
#ifdef INFO_OUTPUT
 printf("tau: %.4f \n",p[i].tau);
#endif

 p[i].scat=p[i].scat+1; 
 /* p[i].abs=0; p[i].scat=1;*/

 return;
}

/***** function fix_temp  ********************************************/

void fix_temp(c)

 struct cell_type c[];
 
{
  int s,j,k;
  int cell;

  /* 
    for(k=0;k<(NLCELLS-1)/2;k++)
     for(j=0;j<NMCELLS-1;j++)
       for(s=0;s<NRCELLS;s++)
       {
         / cell=(int)(s+(k*NMCELLS+j)*NRCELLS);
         cell2=(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS);
  	 printf("+++ %d %d   +++ %d %d %d \n", cell,cell2,s,k,j);     /  
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].nabs=
	   (c[(int)(s+(k*NMCELLS+j)*NRCELLS)].nabs+c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].nabs)/2.;
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].abslum=
	   (c[(int)(s+(k*NMCELLS+j)*NRCELLS)].abslum+c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].abslum)/2.;
       }
  
 
    for(k=((NLCELLS-1)/2);k<NLCELLS-1;k++)
     for(j=0;j<NMCELLS-1;j++)
       for(s=0;s<NRCELLS;s++)
       {
	 /	 cell=(int)(s+(k*NMCELLS+j)*NRCELLS); 
	    cell2=(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS);
	    printf("--- %d %d   --- %d %d %d \n", cell,cell2,s,k,j); /        
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].nabs=c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].nabs;
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].abslum=c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].abslum;
       }
  
  for(s=0;s<NRCELLS;s++)
    for(k=0;k<NLCELLS;k++)
      for(j=0;j<NMCELLS;j++)
	{
	  cell=(int)(s+(k*NMCELLS+j)*NRCELLS);
	  if (c[cell].mass>0) c[cell].temp=solve_ree(c[cell].abslum,c[cell].mass,0,1111111,mean_op_ptr);
	}
  

  */

    for(k=0;k<(NLCELLS-1)/2;k++)
     for(j=0;j<NMCELLS-1;j++)
       for(s=0;s<NRCELLS;s++)
       {
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].temp=
	   (c[(int)(s+(k*NMCELLS+j)*NRCELLS)].temp+c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].temp)/2.;
       }
  
    for(k=((NLCELLS-1)/2);k<NLCELLS-1;k++)
     for(j=0;j<NMCELLS-1;j++)
       for(s=0;s<NRCELLS;s++)
       {	   
	 c[(int)(s+(k*NMCELLS+j)*NRCELLS)].temp=c[(int)(s+((NLCELLS-k-2)*NMCELLS+j)*NRCELLS)].temp;
       }
  
 

#ifdef SPH_TREE


  printf("fix temp not good for SpH... solve ree\n"); exit(0);

#endif
  return;
}

 

