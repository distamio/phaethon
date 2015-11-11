
/*

Functions to find photon frequency and direction depending on the radiation field 

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"
#include "params.h"
#include "phaethon.params.h"
#include "constants.h"
#include "functions.h"
#include "inout.h"
#include "photons.h"
#include "params.grid.h"

/********************* mphotons.c ************************/

#define X1       0.01
#define X2       30

/******** assign_photon_params3 *************************/

/* uniform external radiation field 3D case - non-spherical symmetry      
   photons ejected from the surface of a sphere of radius R_MAX          
   this routine calculates initial position and direction of the photons 
   also initializes other parameters */

void init_isrf_photons_3D(i,p,star,s)

   long i;
   struct photon_packet_type p[];
   struct star_type star[];
   int s;

{

double R1,R2,R3,R4,R5;
double costheta,costheta2,phi,phi2,d1x,d1y,d1z,d2x,d2y,d2z,sintheta,sintheta2;

   R1=(double) (rand())/RAND_MAX;
   R2=(double) (rand())/RAND_MAX;
   R3=(double) (rand())/RAND_MAX;
   R4=(double) (rand())/RAND_MAX;
   R5=(double) (rand())/RAND_MAX;
 
   costheta=1-2*R1;

   phi=2*M_PI*R2;
   sintheta=sin(acos(costheta)); 

   p[i].x=-star[s].radius*sintheta*cos(phi);
   p[i].y=-star[s].radius*sintheta*sin(phi);
   p[i].z=-star[s].radius*costheta;
   p[i].r= star[s].radius;

   costheta2=sqrt(R3); 

   sintheta2=sin(acos(costheta2)); 
   phi2=2*M_PI*R4;

   d1x=sintheta2*cos(phi2);
   d1y=sintheta2*sin(phi2);
   d1z=costheta2;

   d2x=d1x;
   d2y=costheta*d1y+sintheta*d1z;
   d2z=-sintheta*d1y+costheta*d1z;

   p[i].kx=sin(phi)*d2x+cos(phi)*d2y;
   p[i].ky=-cos(phi)*d2x+sin(phi)*d2y;
   p[i].kz=d2z;

   p[i].tau=-log(R5);    
   p[i].abs=0;
   p[i].scat=0;
   p[i].id=1;
   p[i].flag=0;
   p[i].weight=1.;
   p[i].multiplicity=1.;

return;
}



/******** assign_photon_params2 ******************/

/* uniform radiation 1D case - spherical symmetry*/

void init_isrf_photons_1D(i,p,star,s)

   long i;
   struct photon_packet_type p[];
   struct star_type star[];
   int s;

{

double R1,R2;
double aux;

   R1=(double) (rand())/RAND_MAX;
   R2=(double) (rand())/RAND_MAX;
  
   p[i].x=0;
   p[i].y=0;
   p[i].z=star[s].radius;
   p[i].r=star[s].radius;
   p[i].tau=-log(R1);
   
   p[i].kz=-sqrt(R2);    

   aux=sin(acos(p[i].kz));
   p[i].kx=aux;
   p[i].ky=0;
      
   p[i].abs=0;
   p[i].scat=0;
   p[i].id=1;
   p[i].flag=0;
   p[i].weight=1.;
   p[i].multiplicity=1.;
 

return;
}



/******** init_photon_params (initialize photon packets) ***************/

/* point star at (STAR_X,STAR_Y,STAR_Z) */

void init_star_photons(i,p,star,s)

     struct star_type star[];
     int s; 
     long i;
     struct photon_packet_type p[];
        
{
  
  double R1,R2,R3;
  double aux;
  
   R1=(double) (rand())/RAND_MAX;
   R2=(double) (rand())/RAND_MAX;
   R3=(double) (rand())/RAND_MAX;

 
   p[i].kz=2.*R2-1; 
   aux=sin(acos(p[i].kz));
   p[i].kx=aux*cos(2.*M_PI*R3); 
   p[i].ky=aux*sin(2.*M_PI*R3); 
  
   if (star[s].type==0) /* start from R_DUST */
   {
   p[i].x=star[s].x+p[i].kx*star[s].r_dust; 
   p[i].y=star[s].y+p[i].ky*star[s].r_dust;
   p[i].z=star[s].z+p[i].kz*star[s].r_dust;
	}
 
    if (star[s].type==3) /* start from star */
    {
   p[i].x=star[s].x;
   p[i].y=star[s].y;
   p[i].z=star[s].z; 
	}

   p[i].r=pow(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z,0.5);   
    if (star[s].type==0)  p[i].dist_to_star=star[s].r_dust;
	 if (star[s].type==3)  p[i].dist_to_star=0;
   p[i].tau=-log(R1);
   p[i].abs=0;
   p[i].scat=0;
   p[i].id=1;
   p[i].flag=0;
   p[i].weight=1.;
   p[i].multiplicity=1.;
   /* printf("x %.2e y %.2e  z %.2e kx %.2e  ky %.2e  kx %.2e  tau %.2e\n", p[i].x, p[i].y, p[i].z, p[i].kx, p[i].ky, p[i].kz, p[i].tau); */ 
return;
}


/******** Function calculate frequency *******************************/

/* it's the same in x-space for every star                           */

void find_xtable(xtable_ptr,star,s)

    struct xtable_type *xtable_ptr;
    struct star_type star[];
    int s;
  
{
  int i;
  double df;
  double temp=0.;    /* auxiliary table */
  double dlum=0.;          /* luminosity of each packet */
  double bb_constant;
  char filename[20];
  double npackets=0.;
  int fcounter=0;
  
  FILE *fp;
  
  /* i don't need to know the luminosity to distribute the photons */
  
  bb_constant=1.;
  
  
/* create a table with the frequencies (fpoints of them) I want */
/* or, alternatively, I can read a  table with the frequencies  */
 

/* set the frequency space */
/*
df=(X2-X1)/FPOINTS; 

f=df;

for (i=0;i<FPOINTS;i++)
    { 
    (*xtable_ptr).freq[i]=f;
    f=f+df;
    (*xtable_ptr).num[i]=0;
    temp[i]=0;
    } 
*/
 
 
/* set the frequency space [logarithmic intervals]
 */
 df=log10(X2/X1)/FPOINTS; 

 (*xtable_ptr).freq[0]=X1;
 
 for (i=1;i<FPOINTS;i++)
   { 
     (*xtable_ptr).freq[i]=pow(10,(log10((*xtable_ptr).freq[i-1])+df));     
     (*xtable_ptr).num[i]=0;   
   } 

 
 /* calculate photon packet luminosity */
 
 dlum=1./star[s].npackets;
 
 /* calculate number of photon packets per frequency */
 
 /* freq[0]= v1 */
 
 (*xtable_ptr).num[0]=
   (long int) ((simp_b(0.,(((*xtable_ptr).freq[0]+(*xtable_ptr).freq[1])/2.),
		       bb_constant))/dlum);
 
 temp=(*xtable_ptr).num[0]; 

 for (i=1;i<=FPOINTS-2;i++)
   { 
     
     (*xtable_ptr).num[i]=
       (long int) (((simp_b(0.,(((*xtable_ptr).freq[i]+(*xtable_ptr).freq[i+1])/2.),bb_constant))/dlum)-temp);
     temp=temp+(*xtable_ptr).num[i]; 
   }


 (*xtable_ptr).num[FPOINTS-1]=0; /* (long int) (NPACKETS-temp) ; */ /* not important */
 
 
 
 for (i=0;i<=FPOINTS-1;i++)
   
   {
     (*xtable_ptr).multiplicity[i]=1; 
   }
 
 /* print number of photon packets at each frequency */

 sprintf(filename,"fnum.star.%d.dat",s);

 fp = fopen(filename, "w");
 
 fprintf(fp, "fnum.dat : Freq / number of photon packets/ fpoint \n");
 
 for (i=0;i<=FPOINTS-1;i++)
   fprintf(fp,"%.6f %ld %5d \n", (*xtable_ptr).freq[i],(*xtable_ptr).num[i],i );
 
 fclose(fp);
 
 
 /* convert to cgs */
 for (i=0;i<FPOINTS;i++)
   
   (*xtable_ptr).freq[i]=(k_CONST/h_CONST)*(*xtable_ptr).freq[i]*star[s].temp;
 

 sprintf(filename,"fnum.star.cgs%d.dat",s);


 fp = fopen(filename, "w");
 
 fprintf(fp, "fnum.cgs.dat : Freq (hz) / number of photon packets/ lambda (micron) / fpoint / tot photons so far\n");
 

 for (i=0;i<=FPOINTS-1;i++)
   { 
     npackets=npackets+(*xtable_ptr).num[i];
     fprintf(fp,"%.4e %10ld %.4e %5d %12.0f \n", 
    	     (*xtable_ptr).freq[i],(*xtable_ptr).num[i],1e4*c_CONST/(*xtable_ptr).freq[i], i,npackets);
   } 
 fprintf(fp,"\n(source %d) Total Number of Photons %.0f\n",s, npackets);
 
 
 fclose(fp);
 
 
 return;
 
}



/******** Function calculate frequency and luminosity *******************************/

void assign_photon_freq_lum(i,k_ptr,p,xtable_ptr,star,s,pk_fpoint_ptr,pk_num_ptr)

   long i;
   int *k_ptr;
   struct photon_packet_type p[];
   struct xtable_type *xtable_ptr;
   struct star_type star[];
   int s;
   int *pk_fpoint_ptr;
   long *pk_num_ptr;
{


/* assign each photon packet a frequency */

     while  ((*xtable_ptr).num[(*k_ptr)]==0) {(*k_ptr)=(*k_ptr)+1;}  

     if ((*xtable_ptr).num[(*k_ptr)]>0)

       {       
	 p[i].f=(*xtable_ptr).freq[(*k_ptr)];
         p[i].lum=star[s].lum/(star[s].npackets*(*xtable_ptr).multiplicity[(*k_ptr)]);
         p[i].starlum=star[s].lum/SOLAR_LUM;

	 (*xtable_ptr).num[(*k_ptr)]=(*xtable_ptr).num[(*k_ptr)]-1;
	 *pk_fpoint_ptr=(*k_ptr);
	 *pk_num_ptr=(*xtable_ptr).num[(*k_ptr)]; 

        } 

return;

}

/******** Function do_diagnostics *******************************/

void do_diagnostics(p,star,s)

   struct photon_packet_type p[];
   struct star_type star[];
   int s;
  
{

int i;
int k=0;

double tau_mean=0., kx_mean=0., ky_mean=0., kz_mean=0.;
double d_tau_mean=0., d_kx_mean=0., d_ky_mean=0., d_kz_mean=0.;

for (i=0;i<star[s].npackets;i++)
  {

   kx_mean=kx_mean+p[i].kx;
   ky_mean=ky_mean+p[i].ky;
   kz_mean=kz_mean+p[i].kz;
   tau_mean=tau_mean+p[i].tau;

   if ((p[i].kx*p[i].kx+p[i].ky*p[i].ky+p[i].kz*p[i].kz) ==1) k=1;
 }

   kx_mean=kx_mean/star[s].npackets;
   ky_mean=ky_mean/star[s].npackets;
   kz_mean=kz_mean/star[s].npackets;
   tau_mean=tau_mean/star[s].npackets;

   for (i=0;i<star[s].npackets;i++)
  {
   d_kx_mean=d_kx_mean+pow(kx_mean-p[i].kx,2);
   d_ky_mean=d_ky_mean+pow(ky_mean-p[i].ky,2);
   d_kz_mean=d_kz_mean+pow(kz_mean-p[i].kz,2);
   d_tau_mean=d_tau_mean+pow(tau_mean-p[i].tau,2);
  }

  d_kx_mean=sqrt(d_kx_mean/pow(star[s].npackets,2));
  d_ky_mean=sqrt(d_ky_mean/pow(star[s].npackets,2));
  d_kz_mean=sqrt(d_kz_mean/pow(star[s].npackets,2));
  d_tau_mean=sqrt(d_tau_mean/pow(star[s].npackets,2));

   printf("\nDoing diagnostics...\n");
   printf("<kx> =% .5f +- % .5f \n", kx_mean,d_kx_mean);
   printf("<ky> =% .5f +- % .5f \n", ky_mean,d_ky_mean);
   printf("<kz> =% .5f +- % .5f \n", kz_mean,d_kz_mean);
   printf("<tau>=% .5f +- % .5f \n",tau_mean,d_tau_mean);
   printf("Is |k|=1 (1=TRUE)? %d\n\n",k);

return;

}


/******** Function calculate frequency ISRF *******************************/

/* it's the same in x-space for every star                           */

void find_xtable_isrf(xtable_ptr,bisrf_ptr,star,s)

    struct xtable_type *xtable_ptr;  
    struct btype *bisrf_ptr;
    struct star_type star[];
    int s;

{

  int i;
  double df,temp=0;
  double dlum;          /* luminosity of each packet */
  double f_initial;
  double f_final;
  double i_nu_dnu=0.;
  double f1,f2;
  double Y1=(c_CONST/((*bisrf_ptr).lambda[(*bisrf_ptr).dim-1]*1.e-4));    /* frequency1 */
  double Y2=(c_CONST/((*bisrf_ptr).lambda[0]*1.e-4));                      /* frequency2 */
  double integral=0.;
  double npackets=0;
  char filename[20];
  int fcounter=0;
FILE *fp;
 
/* UNITS: I(nu) dnu =  0.00432656 erg s-1 cm-2 black */
/* 4.91928e-03 endre */

f_initial=c_CONST/((*bisrf_ptr).lambda[(*bisrf_ptr).dim-1]*1.e-4);
f_final=c_CONST/((*bisrf_ptr).lambda[0]*1.e-4);

(*bisrf_ptr).int_total=simp_isrf(f_initial,f_final,bisrf_ptr,24);
printf("\n::Intensity integrated over nu     :%.10e\n",(*bisrf_ptr).int_total);

/*
i1=simp_isrf(c_CONST/(bisrf.lambda[(*bisrf_ptr).dim-1]*1.e-4),c_CONST/(400*1.e-4),bisrf,20);
i2=simp_isrf(c_CONST/(400*1.e-4),c_CONST/(5*1.e-4),bisrf,20);
i3=simp_isrf(c_CONST/(5*1.e-4),c_CONST/(bisrf.lambda[0]*1.e-4),bisrf,20);
*/

/* set the frequency space */
/* 
df=(Y2-Y1)/FPOINTS;  
f=df; y
for (i=0;i<FPOINTS;i++) 
    {  
     (*xtable_ptr).freq[i]=f; 
     f=f+df; 
     (*xtable_ptr).num[i]=0; 
     temp[i]=0; 
 } 
*/

/* set the frequency space [logarithmic intervals] */


df=log10(Y2/Y1)/FPOINTS; 

(*xtable_ptr).freq[0]=Y1;

for (i=1;i<FPOINTS;i++)
    { 

     (*xtable_ptr).freq[i]=pow(10,(log10((*xtable_ptr).freq[i-1])+df));     
     (*xtable_ptr).num[i]=0;
    
    }    

/* calculate total luminosity (again) */

f2=((*xtable_ptr).freq[0]+(*xtable_ptr).freq[1])/2.;
f2=(*xtable_ptr).freq[1];//

integral=simp_isrf(Y1,f2,bisrf_ptr,9);

for (i=1;i<=FPOINTS-2;i++)
    {
    f1=f2; 
    f2=( (*xtable_ptr).freq[i]+(*xtable_ptr).freq[i+1])/2.;
     f2=(*xtable_ptr).freq[i+1]; // 
    integral=simp_isrf(f1,f2,bisrf_ptr,9)+integral;
    }

 (*bisrf_ptr).int_total=integral;
 printf("::Intensity integrated over nu (2) :%.10e\n",(*bisrf_ptr).int_total);  
 integral=0.;
 
/* calculate photon packet luminosity  (I(nu) nu/ npackets) */

dlum=(*bisrf_ptr).int_total/star[s].npackets; 

/* calculate number of photon packets per frequency */

f2=((*xtable_ptr).freq[0]+(*xtable_ptr).freq[1])/2.;

f2=(*xtable_ptr).freq[1]; //
integral=simp_isrf(Y1,f2,bisrf_ptr,9);

(*xtable_ptr).num[0]=(long int) (integral/dlum);

 npackets=npackets+(*xtable_ptr).num[0];

temp=(*xtable_ptr).num[0];


for (i=1;i<=FPOINTS-2;i++)
    {

    f1=f2;
    f2=( (*xtable_ptr).freq[i]+(*xtable_ptr).freq[i+1])/2.;
    f2=(*xtable_ptr).freq[i+1]; // 
    

    integral=simp_isrf(f1,f2,bisrf_ptr,9)+integral;

    (*xtable_ptr).num[i]=(long int)((integral/dlum)-temp);
    temp=temp+(*xtable_ptr).num[i];
    npackets=npackets+(*xtable_ptr).num[i];
    if (npackets>star[s].npackets) { fcounter=fcounter+(*xtable_ptr).num[i];(*xtable_ptr).num[i]=0;}
   
    }



if  (fcounter>100) {printf("ERROR (1) photons.c :: fcounter=%.d\n",fcounter);exit(0);} 
 
 (*xtable_ptr).num[FPOINTS-1]= (long int) (star[s].npackets-temp) ;
 
 if ( (*xtable_ptr).num[FPOINTS-1]<0 ) (*xtable_ptr).num[FPOINTS-1]=0;
 
 if ( (*xtable_ptr).num[FPOINTS-1]<-10 ) {printf("ERROR (2) photons.c\n");exit(0);}
 
 npackets=npackets+(*xtable_ptr).num[FPOINTS-1];
 
 i_nu_dnu=integral;
 printf("::Intensity integrated over nu (3) :%.10e\n",i_nu_dnu);

/* arrange for more photons to be emitted at a specific frequency */
 
find_xtable_multiplicity(xtable_ptr) ;

 for (i=0;i<FPOINTS;i++)
    
 (*xtable_ptr).num[i]= (*xtable_ptr).num[i]*(*xtable_ptr).multiplicity[i];


/* print number of photon packets at each frequency */

sprintf(filename,"fnum.star.%d.dat",s);
 
fp = fopen(filename, "w");

fprintf(fp, "fnum.dat : Freq / number of photon packets\n");

for (i=0;i<=FPOINTS-1;i++)
fprintf(fp,"%.4e %ld\n", (*xtable_ptr).freq[i],(*xtable_ptr).num[i]);


fclose(fp);


 npackets=0;

sprintf(filename,"fnum.star.cgs.%d.dat",s);
 
 fp = fopen(filename, "w");
 
 
fprintf(fp, "fnum.cgs.dat : Freq (hz) / number of photon packets/ lambda (micron) / fpoint / tot photons so far\n");
 
 for (i=0;i<=FPOINTS-1;i++)
   { 
     npackets=npackets+(*xtable_ptr).num[i];
     fprintf(fp,"%.4e %10ld %.4e %5d %12.0f \n", 
 	     (*xtable_ptr).freq[i],(*xtable_ptr).num[i],1e4*c_CONST/(*xtable_ptr).freq[i], i,npackets);
   } 
 fprintf(fp,"\n(source %d) Total Number of Photons %.0f\n",s, npackets);

fclose(fp);


return;

}



/************************************************************************/

void find_xtable_multiplicity(xtable_ptr)

    struct xtable_type *xtable_ptr;  

{
  int i;

  if (MULTIPLE_PHOTONS==1) printf("\nWarning :: Multiplicity HAS NOT BEEN TESTED IN NEW PHAETHON!!!\n");
 
  (*xtable_ptr).multiplicity[0]=1;

 if (MULTIPLE_PHOTONS==1)

   {
     for (i=1;i<=FPOINTS-1;i++)
      {
	
	if (1.e4*c_CONST/(*xtable_ptr).freq[i]<8) 
	  
           (*xtable_ptr).multiplicity[i]=300;
	
	
       /*
           if (1.e4*c_CONST/(*xtable_ptr).freq[i]>1  && 1.e4*c_CONST/(*xtable_ptr).freq[i]<20)
 
           (*xtable_ptr).multiplicity[i]=100;
       */


      else (*xtable_ptr).multiplicity[i]=0;
 
       }
  }

    else
      {

       for (i=0;i<=FPOINTS-1;i++)

         {
          (*xtable_ptr).multiplicity[i]=1; 
         }
   
      } 

  return;

}


/************ assing opacity for a photon *************************/

void assign_photon_opacity(p,i,opa_ptr, fpoint)
  
     struct photon_packet_type p[];
     int i;
     struct opacity_table * opa_ptr;
     int fpoint;
{
  p[i].k_abs=(*opa_ptr).abs[fpoint];
  p[i].s_scat=(*opa_ptr).scat[fpoint];
  p[i].albedo=(*opa_ptr).albedo[fpoint];
  
  return;
  
}


/************** assing opacity for a photon *************************/
#ifdef DUST_MIX
void assign_photon_opacity2(p,i,opa_ptr, fpoint)
  
     struct photon_packet_type p[];
     int i;
    struct opacity_table * opa_ptr;
     int fpoint;
{
  p[i].k_abs2=(*opa_ptr).abs[fpoint];
  p[i].s_scat2=(*opa_ptr).scat[fpoint];
  p[i].albedo2=(*opa_ptr).albedo[fpoint];
  
  return;
}  
#endif



















/******** assign_photon_params_box *************************/

/* uniform external radiation field 3D case - non-spherical symmetry      
   photons ejected from the surface of a sphere of radius R_MAX          
   this routine calculates initial position and direction of the photons 
   also initializes other parameters */

void init_photon_params_box(i,p,rsize)

   long i;
   struct photon_packet_type p[];
   double rsize;

{

double R1,R2,R3,R4,R5,R6;
double costheta,costheta2,phi,phi2,sintheta,sintheta2,kx,ky,kz;


   R1=(double) (rand())/RAND_MAX;
   R2=(double) (rand())/RAND_MAX;
   R3=(double) (rand())/RAND_MAX;
   R4=(double) (rand())/RAND_MAX;
   R5=(double) (rand())/RAND_MAX;
   R6=(double) (rand())/RAND_MAX;

  
   /* position on the box:  R1, R2, R6 */

   R1=6.*R1;
   
   if (R1<=1)         {p[i].x=-0.5*rsize; p[i].y=(1-2.*R2)*0.5*rsize; p[i].z=(1-2.*R6)*0.5*rsize;}
   if (R1<=2 && R1>1) {p[i].x=+0.5*rsize; p[i].y=(1-2.*R2)*0.5*rsize; p[i].z=(1-2.*R6)*0.5*rsize;} 
   if (R1<=3 && R1>2) {p[i].y=-0.5*rsize; p[i].x=(1-2.*R2)*0.5*rsize; p[i].z=(1-2.*R6)*0.5*rsize;} 
   if (R1<=4 && R1>3) {p[i].y=+0.5*rsize; p[i].x=(1-2.*R2)*0.5*rsize; p[i].z=(1-2.*R6)*0.5*rsize;} 
   if (R1<=5 && R1>4) {p[i].z=-0.5*rsize; p[i].y=(1-2.*R2)*0.5*rsize; p[i].x=(1-2.*R6)*0.5*rsize;} 
   if (R1>5)          {p[i].z=+0.5*rsize; p[i].y=(1-2.*R2)*0.5*rsize; p[i].x=(1-2.*R6)*0.5*rsize;} 
   
   p[i].r=sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);  

   costheta=p[i].z/p[i].r; 

   sintheta=sin(acos(costheta));

   phi=atan2(p[i].y,p[i].x); if  (p[i].phi<0) p[i].phi=p[i].phi+2.*M_PI; 
   
   /*    
   if (R1<=1)         {phi=0;  theta=M_PI/2.;} 
   if (R1<=2 && R1>1) {phi=0;  theta=-M_PI/2.;} 
   if (R1<=3 && R1>2) {phi=0;  theta=M_PI/2.;} 
   if (R1<=4 && R1>3) {phi=0;  theta=-M_PI/2.;} 
   if (R1<=5 && R1>4) {phi=0;  theta=0;} 
   if (R1>5)          {phi=0;  theta=M_PI;} 

   costheta=cos(theta);
   sintheta=pow(1-costheta*costheta,0.5); 
   */
 
   /* direction of photon */
   
   costheta2=sqrt(R3); 
   sintheta2=sin(acos(costheta2));
   phi2=2*M_PI*R4;  
   
   kx=sintheta2*cos(phi2);
   ky=sintheta2*sin(phi2);
   kz=costheta2;  

   /* direction vector coordinates in photons idiosystem */
   /*  
       d1x=sintheta2*cos(phi2);
       d1y=sintheta2*sin(phi2);
       d1z=costheta2;
   */
   /* transform to original system */
   
   /* 
      d2x=d1x;
      d2y=costheta*d1y+sintheta*d1z;
      d2z=-sintheta*d1y+costheta*d1z;
      
      p[i].kx=sin(phi)*d2x+cos(phi)*d2y;
      p[i].ky=-cos(phi)*d2x+sin(phi)*d2y;
      p[i].kz=d2z;
      
   */

   if (R1<=1)          {p[i].kx= kz;  p[i].ky= ky;  p[i].kz=-kx;}
   if (R1<=2 && R1>1)  {p[i].kx=-kz;  p[i].ky=-ky;  p[i].kz= kx;}
   if (R1<=3 && R1>2)  {p[i].kx=-ky;  p[i].ky= kz;  p[i].kz= kx;}
   if (R1<=4 && R1>3)  {p[i].kx =ky;  p[i].ky=-kz;  p[i].kz=-kx;}
   if (R1<=5 && R1>4)  {p[i].kx= kx;  p[i].ky= ky;  p[i].kz= kz;}
   if (R1>5)           {p[i].kx=-kx;  p[i].ky=-ky;  p[i].kz=-kz;}

   if (R1<=1)         {if (p[i].kx<0) printf("er1\n");exit;} 
   if (R1<=2 && R1>1) {if (p[i].kx>0) printf("er2\n");exit;}
   if (R1<=3 && R1>2) {if (p[i].ky<0)  printf("er3\n");exit;}
   if (R1<=4 && R1>3) {if (p[i].ky>0)  printf("er4\n");exit;}
   if (R1<=5 && R1>4) {if (p[i].kz<0)  printf("er5\n");exit;}
   if (R1>5)          {if (p[i].kz>0)  printf("er6\n");exit;}


   p[i].tau=-log(R5);    
   p[i].abs=0;
   p[i].scat=0;
   p[i].id=1;
   p[i].flag=0;
   p[i].weight=1.; 
   p[i].multiplicity=1;

return;
}













/******** Function calculate frequency STAR SED *******************************/

/* it's the same in x-space for every star                           */

void find_xtable_star_sed(xtable_ptr,star,s)

    struct xtable_type *xtable_ptr;  
    struct star_type star[];
    int s;

{

  int i;
  double df,temp=0;
  double dlum;          /* luminosity of each packet */
  double f_initial;
  double f_final;
  double i_nu_dnu=0.;
  double f1,f2;
  double Y1=(c_CONST/(star[s].lambda[star[s].sed_dim-1]*1.e-4));    /* frequency1 */
  double Y2=(c_CONST/(star[s].lambda[0]*1.e-4));                      /* frequency2 */
  double integral=0.;
  double npackets=0;
  char filename[20];
  int fcounter=0;
FILE *fp;
 
/* UNITS: I(nu) dnu =  0.00432656 erg s-1 cm-2 black */
/* 4.91928e-03 Andre */

f_initial=c_CONST/(star[s].lambda[star[s].sed_dim-1]*1.e-4);
f_final=c_CONST/(star[s].lambda[0]*1.e-4);

star[s].int_total=simp_star_sed(f_initial,f_final,star,s,24);
printf("Intensity integrated over nu     :%.10e\n",star[s].int_total);


/* set the frequency space */
/* 
df=(Y2-Y1)/FPOINTS;  
f=df; y
for (i=0;i<FPOINTS;i++) 
    {  
     (*xtable_ptr).freq[i]=f; 
     f=f+df; 
     (*xtable_ptr).num[i]=0; 
     temp[i]=0; 
 } 
*/

/* set the frequency space [logarithmic intervals] */


df=log10(Y2/Y1)/FPOINTS; 

(*xtable_ptr).freq[0]=Y1;

for (i=1;i<FPOINTS;i++)
    { 

     (*xtable_ptr).freq[i]=pow(10,(log10((*xtable_ptr).freq[i-1])+df));     
     (*xtable_ptr).num[i]=0;
    
    }    

/* calculate total luminosity (again) */

f2=((*xtable_ptr).freq[0]+(*xtable_ptr).freq[1])/2.;

 integral=simp_star_sed(Y1,f2,star,s,9); 

for (i=1;i<=FPOINTS-2;i++)
    {
    f1=f2; 
    f2=( (*xtable_ptr).freq[i]+(*xtable_ptr).freq[i+1])/2.;
      integral=simp_star_sed(f1,f2,star,s,9)+integral;
    }

 star[s].int_total=integral;
 printf("Intensity integrated over nu (2) :%.10e\n",star[s].int_total);  
 integral=0.;
 
/* calculate photon packet luminosity  (I(nu) nu/ npackets) */

dlum=star[s].int_total/star[s].npackets; 

/* calculate number of photon packets per frequency */

f2=((*xtable_ptr).freq[0]+(*xtable_ptr).freq[1])/2.;

 integral=simp_star_sed(Y1,f2,star,s,9); 

(*xtable_ptr).num[0]=(long int) (integral/dlum);

 npackets=npackets+(*xtable_ptr).num[0];

temp=(*xtable_ptr).num[0];


for (i=1;i<=FPOINTS-2;i++)
    {

    f1=f2;
    f2=( (*xtable_ptr).freq[i]+(*xtable_ptr).freq[i+1])/2.;
    
    integral=simp_star_sed(f1,f2,star,s,9)+integral; 

    (*xtable_ptr).num[i]=(long int)((integral/dlum)-temp);
    temp=temp+(*xtable_ptr).num[i];
    npackets=npackets+(*xtable_ptr).num[i];
    if (npackets>star[s].npackets) { fcounter=fcounter+(*xtable_ptr).num[i];(*xtable_ptr).num[i]=0;}
   
    }



if  (fcounter>100) {printf("ERROR (1) photons.c :: fcounter=%.d\n",fcounter);exit(0);} 
 
 (*xtable_ptr).num[FPOINTS-1]= (long int) (star[s].npackets-temp) ;
 
 if ( (*xtable_ptr).num[FPOINTS-1]<0 ) (*xtable_ptr).num[FPOINTS-1]=0;
 
 if ( (*xtable_ptr).num[FPOINTS-1]<-10 ) {printf("ERROR (2) photons.c\n");exit(0);}
 
 npackets=npackets+(*xtable_ptr).num[FPOINTS-1];
 
 i_nu_dnu=integral;
 printf("Intensity integrated over nu (3) :%.10e\n",i_nu_dnu);

/* arrange for more photons to be emitted at a specific frequency */
 
find_xtable_multiplicity(xtable_ptr) ;

 for (i=0;i<FPOINTS;i++)
    
 (*xtable_ptr).num[i]= (*xtable_ptr).num[i]*(*xtable_ptr).multiplicity[i];


/* print number of photon packets at each frequency */

sprintf(filename,"fnum.star.%d.dat",s);
 
fp = fopen(filename, "w");

fprintf(fp, "fnum.dat : Freq / number of photon packets\n");

for (i=0;i<=FPOINTS-1;i++)
fprintf(fp,"%.4e %ld\n", (*xtable_ptr).freq[i],(*xtable_ptr).num[i]);


fclose(fp);


 npackets=0;

sprintf(filename,"fnum.star.cgs.%d.dat",s);
 
 fp = fopen(filename, "w");
 
 
fprintf(fp, "fnum.cgs.dat : Freq (hz) / number of photon packets/ lambda (micron) / fpoint / tot photons so far\n");
 
 for (i=0;i<=FPOINTS-1;i++)
   { 
     npackets=npackets+(*xtable_ptr).num[i];
     fprintf(fp,"%.4e %10ld %.4e %5d %12.0f \n", 
 	     (*xtable_ptr).freq[i],(*xtable_ptr).num[i],1e4*c_CONST/(*xtable_ptr).freq[i], i,npackets);
   } 
 fprintf(fp,"\n(source %d) Total Number of Photons %.0f\n",s, npackets);

fclose(fp);


return;

}

