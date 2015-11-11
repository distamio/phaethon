/*
   Functions that Calculate the SED of the system */

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




/**************** function create_fbins ****************************/

/* creates logarithmic output frequency bins */

void create_fbins(fbin_ptr,opacity)

  struct fbin_type *fbin_ptr;
  struct opacity_type opacity;
{

int k;
double dlambda;


(*fbin_ptr).lambda[0]=SL1; /* starting lambda (microns) */
dlambda=log10(SL2/SL1)/NFBINS; /* logarithic fbins      */

 for (k=1;k<=NFBINS;k++)
  {
   (*fbin_ptr).lambda[k]=pow(10,(log10((*fbin_ptr).lambda[k-1])+dlambda));

  }
 
 for (k=0;k<=NFBINS;k++)
   {
     (*fbin_ptr).lambda1[k]=pow(10,-dlambda/2.)*(*fbin_ptr).lambda[k];
     (*fbin_ptr).lambda2[k]=pow(10,dlambda/2.)*(*fbin_ptr).lambda[k];
   }

 init_fbins(fbin_ptr,opacity);

 return;
}

/**************** function init_fbins ****************************/

/* initialiaze frequency bins */
void init_fbins(fbin_ptr,opacity)

  struct fbin_type *fbin_ptr; struct opacity_type opacity;
{
  
  int k,id,l,m;
  
  for (k=0;k<=NFBINS;k++)
    {
      (*fbin_ptr).abs[k]=k_abs(1e4*c_CONST/(*fbin_ptr).lambda[k],opacity);
      (*fbin_ptr).scat[k]=sigma_scat(1e4*c_CONST/(*fbin_ptr).lambda[k],opacity);
      (*fbin_ptr).opacity[k]=k_abs(1e4*c_CONST/(*fbin_ptr).lambda[k],opacity)+
	sigma_scat(1e4*c_CONST/(*fbin_ptr).lambda[k],opacity);
      
      /* printf ("%.4f %.4f\n",(*fbin_ptr).lambda[k],(*fbin_ptr).opacity[k]); */
      
      for (id=0;id<=MAX_ID;id++)
	{
	  for (l=0;l<NLBINS;l++)
	    for (m=0;m<NMBINS;m++)
	      {
		(*fbin_ptr).lum[k][id][l][m]=0.;
		(*fbin_ptr).npackets[k][id][l][m]=0.;
	      }
	}
    }
  return;
 
}

/**************** function 
 if (p[i].id==4) printf("%d\n",p[i].id);
create_dbins ****************************/

/* define direction bins */

void create_dbins(dbin_ptr)
     
     struct dbin_type *dbin_ptr;
{
  
  /* theta and phi must be transformed in RAD in this function */

  int l,m;

   for (l=0;l<NLBINS;l++) 
     for (m=0;m<NMBINS;m++) 
       {
	 switch (m)
	   {
	   case 0:(*dbin_ptr).phi[l][m]=PHI0;break;
	   case 1:(*dbin_ptr).phi[l][m]=PHI1;break;
	   case 2:(*dbin_ptr).phi[l][m]=PHI2;break;
	   case 3:(*dbin_ptr).phi[l][m]=PHI3;break;
	   case 4:(*dbin_ptr).phi[l][m]=PHI4;break;
	   case 5:(*dbin_ptr).phi[l][m]=PHI5;break;
	   case 6:(*dbin_ptr).phi[l][m]=PHI6;break;
	   case 7:(*dbin_ptr).phi[l][m]=PHI7;break;
	   case 8:(*dbin_ptr).phi[l][m]=PHI8;break;
	   case 9:(*dbin_ptr).phi[l][m]=PHI9;break;
	   }
     
       }
  

   for (l=0;l<NLBINS;l++) 
     for (m=0;m<NMBINS;m++) 
       {	 
	 switch (l)
	   {
	   case 0:(*dbin_ptr).theta[l][m]=THETA0;break;
	   case 1:(*dbin_ptr).theta[l][m]=THETA1;break;
	   case 2:(*dbin_ptr).theta[l][m]=THETA2;break;
	   case 3:(*dbin_ptr).theta[l][m]=THETA3;break;
	   case 4:(*dbin_ptr).theta[l][m]=THETA4;break;
	   case 5:(*dbin_ptr).theta[l][m]=THETA5;break;
	   case 6:(*dbin_ptr).theta[l][m]=THETA6;break;
	   case 7:(*dbin_ptr).theta[l][m]=THETA7;break;
	   case 8:(*dbin_ptr).theta[l][m]=THETA8;break;
	   case 9:(*dbin_ptr).theta[l][m]=THETA9;break;
	   }
       }

   printf("\n l m  theta   phi\n");
  


   for (l=0;l<NLBINS;l++)
     for (m=0;m<NMBINS;m++)
       
       printf(" %d %d %6.1f %5.1f\n",l,m, (*dbin_ptr).theta[l][m],(*dbin_ptr).phi[l][m] );
   

   for (l=0;l<NLBINS;l++)
     for (m=0;m<NMBINS;m++)
       {
	 (*dbin_ptr).theta[l][m]=(*dbin_ptr).theta[l][m]*M_PI/180.;
	 (*dbin_ptr).phi[l][m]=(*dbin_ptr).phi[l][m]*M_PI/180.;
	   }   
   

   calculate_domega(dbin_ptr);

   return;
   
}


/*************** function place_into_fbins ***********************/
void place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s)

     struct photon_packet_type p[];
     int i;
     struct fbin_type *fbin_ptr;
     double *total_fbin_ptr;
     int l, m;
     struct star_type star[];
     int s;
     
{
  int k=0;
  double lambda;       


  if  (star[s].type==3 && p[i].r==star[s].z)
    
    return;
      
  else
    {
      /* dlambda=log10(SL2/SL1)/NFBINS;*/ /* logarithic fbins      */
      
      lambda=(c_CONST/p[i].f)*1.e4;
      
      k=binary_search_fbin (lambda,0, NFBINS,fbin_ptr);
 
      
      if (k>=0)     
	{
	  p[i].n=k;
	  
	  (*fbin_ptr).lum[p[i].n][0][l][m]=(*fbin_ptr).lum[p[i].n][0][l][m]+p[i].weight*p[i].lum;
	  (*fbin_ptr).npackets[p[i].n][0][l][m]=(*fbin_ptr).npackets[p[i].n][0][l][m]+p[i].weight;
	 
     if (!(p[i].id==3 && (star[s].type==1 || star[s].type==2) && (ISRF_IGNORE_SCATTER==1)))
 
	  if (p[i].origin==1)
	    {
	      (*fbin_ptr).lum[p[i].n][2][l][m]=(*fbin_ptr).lum[p[i].n][2][l][m]+p[i].weight*p[i].lum;
	      (*fbin_ptr).npackets[p[i].n][2][l][m]=(*fbin_ptr).npackets[p[i].n][2][l][m]+p[i].weight;
	    }
	  
	  (*fbin_ptr).lum[p[i].n][p[i].id][l][m]=(*fbin_ptr).lum[p[i].n][p[i].id][l][m]+p[i].weight*p[i].lum;	  
	  (*fbin_ptr).npackets[p[i].n][p[i].id][l][m]=(*fbin_ptr).npackets[p[i].n][p[i].id][l][m]+p[i].weight; 
 
	  
	  (*total_fbin_ptr)=(*total_fbin_ptr)+p[i].weight;
	} 
	

    }
  
  return;
  
}


/******************** modified binary search **************************/

int binary_search_fbin (lambda, low, high,fbin_ptr)


double lambda;
int low;
int high;
struct fbin_type *fbin_ptr;

{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (lambda>=(*fbin_ptr).lambda1[middle] && lambda<(*fbin_ptr).lambda2[middle]) 

    return middle;

 else if (lambda<(*fbin_ptr).lambda1[middle])

         high=middle-1;

 else low=middle+1;

 }

return (-1);

}


/*************** function place_into_dbins ***********************/

void place_into_dbins(p,i,dbin_ptr,total_dbin_ptr)

struct photon_packet_type p[];
int i;
struct dbin_type *dbin_ptr;
double *total_dbin_ptr;

{
double phi,cos_theta;
double aux;

      /* find phi */

  if (p[i].kx==0) phi=M_PI/2.; 

 else
   {   aux=atan(p[i].ky/p[i].kx);

         if (p[i].kx<0 && p[i].ky>=0) phi=-aux+M_PI/2.;
        
         else 
              {
		if (p[i].kx<0 && p[i].ky<0) phi=M_PI+aux;

	       else 

               if (p[i].kx>0 && p[i].ky<0) phi=2.*M_PI+aux; 

               else 

		 phi=aux;
     	     }
     }
	 /* find theta */

         cos_theta=(p[i].kz);       



	 /***** NEED TO CHANGE THE ABOVE!!!! ********/

return;

}


/******************** modified binary search **************************/

int binary_search_bbin (r, low, high, bbin_ptr)


double r;
int low;
int high;
struct bbin_type *bbin_ptr;

{

int middle;

 while (low<=high) {

 middle= (low+high)/2;

 if (r>=(*bbin_ptr).b[middle] && r<(*bbin_ptr).b[middle+1]) 

    return middle;

 else if (r<(*bbin_ptr).b[middle])

         high=middle-1;

 else low=middle+1;

 }

printf("search key not found!\n"); 
exit(0);/* search key not found */

return -1;

}
/********** function bbin intensity **************************************/

/* returns intensity for a given r (spherically symmetric case) and fbin*/

#ifdef DO_ISOMAP
void calculate_isogrid_intensity(i,j,bbin_ptr,isogrid_ptr,background_ptr,startype)
     
     double *background_ptr;
     int i,j;
     struct bbin_type *bbin_ptr;
     struct isogrid_type *isogrid_ptr;
	
     int startype;
{
  int index,f,id;
  double r;
  int l=0,m=0;
  /* bbin.b[NBBINS] is always zero! last cell! */
  
  r=sqrt(pow((*isogrid_ptr).x[i][j],2)+pow((*isogrid_ptr).y[i][j],2));
  


  if (r>(*bbin_ptr).b[NBBINS]) 

    { 
      for (id=0;id<=MAX_ID;id++)
	{
	  for (f=0;f<FS;f++)
	    {
	      if(id==3 || id==4 || id==2) (*isogrid_ptr).intensity[i][j][f][id][l][m]=0.;

	      if (startype==1 || startype==2 || startype==4) (*isogrid_ptr).intensity[i][j][f][id][l][m]=*(background_ptr+f); 
	      /*	     if (STAR==1 && EXTERNAL_UNIFORM_RADIATION_3D==1) correct but need to include lambda in the funtion
			     (*isogrid_ptr).intensity[i][j][f][id]=lblack((*(lambda_obs_ptr+f))*1e-4,STAR_TEMP)*
			     1e-8*((*(lambda_obs_ptr+f))*(*(lambda_obs_ptr+f))/c_CONST)*1e17; */
	    }
	}
    }
  
  else 
    {
      index=binary_search_bbin(r,0,NBBINS-1,bbin_ptr);
      
      for (id=0;id<=MAX_ID;id++)
	{
	  for (f=0;f<FS;f++)
	    {
	      
	      if (index==NBBINS-1) 
		
		
		(*isogrid_ptr).intensity[i][j][f][id][l][m]=(*bbin_ptr).intensity[NBBINS-1][f][id];
	      
	      else
		{
		  
		  
		  (*isogrid_ptr).intensity[i][j][f][id][l][m]=((*bbin_ptr).intensity[index+1][f][id]*(r-(*bbin_ptr).b[index])+ 
							       ((*bbin_ptr).intensity[index][f][id]*((*bbin_ptr).b[index+1]-r)) )/
		    ((*bbin_ptr).b[index+1]-(*bbin_ptr).b[index]);
		}
	      
	    }
	}
    }
  
  
  return;
  
}
#endif

/***** function write_isophotal_map2  *****************************************/

/* cgs units/ spherically symmetric case */

void write_isophotal_map2(bbin_ptr,filename)

     struct bbin_type *bbin_ptr;
     char filename[20];
{
  int i,k,id;
  double b;
  char filename1[20];
  FILE *fp;
  
  for (k=0;k<FS;k++)
    {
      
      for(id=0;id<=MAX_ID;id++)
	{
	  
	  sprintf(filename1,"%s.%d.%d.dat",filename,k,id);
	  
	  /* Open  file */
	  if (( fp = fopen(filename1, "w")) == NULL)
	    {
	      fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename1);
	      return;
	    }
	  
	  
	  fprintf(fp,"@b(cm)/I_nu(b)(MJy)(no_spaces) \n");
	  
	  for (i =0; i<=NBBINS-1; i++)
	    {
	      b=((*bbin_ptr).b[i]+(*bbin_ptr).b[i+1])/2.;     
	      fprintf(fp,"% .3e % .3e \n", b, (*bbin_ptr).intensity[i][k][id]);
	    } 
	  
	  fclose(fp);
	}
    }
  
  return;
}


/******* make_isomap *******************************************************/

/* this is for the special case of spherical symmetry
 */

#ifdef DO_ISOMAP
void make_isogrid_spherical_case(bbin_ptr,isogrid_ptr,background_ptr,startype)

        struct bbin_type *bbin_ptr;
	struct isogrid_type *isogrid_ptr;
	double *background_ptr;	
	
        int startype;
{
  int i,j,k;
  
  
  for (i=0;i<NXPIXELS;i++)
    
    for  (j=0;j<NYPIXELS;j++)
      
      {
	(*isogrid_ptr).x[i][j]=XMAX-(i+1./2)*(*isogrid_ptr).dx;
	(*isogrid_ptr).y[i][j]=YMAX-(j+1./2)*(*isogrid_ptr).dy;
	   
	
	
 	calculate_isogrid_intensity(i,j,bbin_ptr,isogrid_ptr,background_ptr,startype);
      }

 return;
}
#endif

/*************** function place_into_bbins ***********************/
#ifdef SPHERICAL
void place_into_bbins(p,i,bbin_ptr,total_bbin_ptr,fbin_obs_ptr,fbin_obs_min_ptr,fbin_obs_max_ptr,r_max,star,s)
     
     struct photon_packet_type p[];
     int i;
     struct bbin_type *bbin_ptr;
     double *total_bbin_ptr;
     int *fbin_obs_ptr; /* obsrveving lambda bin (if -1 then all) */
     int *fbin_obs_min_ptr;
     int *fbin_obs_max_ptr;
     double r_max;
     struct star_type star[];
     int s;

{
  int k,f;
  
  for(f=0;f<FS;f++)
    for (k=0;k<=NBBINS-1;k++)
   {  
     if  (p[i].b>=r_max) p[i].b=r_max;
     if ( p[i].b>=(*bbin_ptr).b[k] && p[i].b<=(*bbin_ptr).b[k+1] && 
	  ( ( p[i].n>=(*(fbin_obs_min_ptr+f)) && p[i].n<=(*(fbin_obs_max_ptr+f)) ) || (*(fbin_obs_ptr+f))<0) ) 
       
       {        
	 (*bbin_ptr).lum[k][f][0]=(*bbin_ptr).lum[k][f][0]+p[i].weight *p[i].lum;
	 (*bbin_ptr).npackets[k][f][0]=(*bbin_ptr).npackets[k][f][0]+p[i].weight;

     if (!(p[i].id==3 && (star[s].type==1 || star[s].type==2) && (ISRF_IGNORE_SCATTER==1)))

	 if (p[i].origin==1)
	   {
	     (*bbin_ptr).lum[k][f][2]=(*bbin_ptr).lum[k][f][2]+p[i].weight *p[i].lum;
	     (*bbin_ptr).npackets[k][f][2]=(*bbin_ptr).npackets[k][f][2]+p[i].weight;
	   }

	 (*bbin_ptr).lum[k][f][p[i].id]=(*bbin_ptr).lum[k][f][p[i].id]+p[i].weight*p[i].lum;	 
	 (*bbin_ptr).npackets[k][f][p[i].id]=(*bbin_ptr).npackets[k][f][p[i].id]+p[i].weight;
	 *(total_bbin_ptr+f)=*(total_bbin_ptr+f)+p[i].weight;
	 
       }
     
   }
 
 return;

}
#endif

/**************** function create_bbins ****************************/
void create_bbins(bbin_ptr,r_max)
  struct bbin_type *bbin_ptr;
  double r_max;
{

int k;
double dr;

 
     (*bbin_ptr).b[0]=0.;
     
     /* logaritmic
	(*bbin_ptr).b[1]=R_DUST_FACTOR*R_MAX;
	
	dr=log10(R_MAX/(*bbin_ptr).b[1])/NBBINS;
	
	for(k=2;k<=NBBINS;k++){(*bbin_ptr).b[k]=pow(10,(log10((*bbin_ptr).b[k-1])+dr));printf("---%.3e\n",(*bbin_ptr).b[k]);}
	
	(*bbin_ptr).b[NBBINS]=R_MAX;
     */
     
     dr=1.*r_max/NBBINS;
     
     for(k=1;k<=NBBINS;k++) (*bbin_ptr).b[k]=k*dr; 
     
     (*bbin_ptr).b[NBBINS]=r_max;
     
     if (MORE==0)
       
       init_bbins(bbin_ptr);
     
     else 
       
       read_bbins(bbin_ptr);
     
return;
 
}

/**************** function init_bbins ****************************/

void init_bbins(bbin_ptr)

  struct bbin_type *bbin_ptr;
{

  int f,k,id;

  for (id=0;id<=MAX_ID;id++)
    
    for (f=0;f<FS;f++)
      
      for (k=0;k<=NBBINS;k++)
	{
	  (*bbin_ptr).lum[k][f][id]=0.;
	  (*bbin_ptr).npackets[k][f][id]=0;
	  (*bbin_ptr).intensity[k][f][id]=0.;
	}
  
  return;
  
}

/**************** function init_bbins ****************************/

void null_bbins_npackets(bbin_ptr)

  struct bbin_type *bbin_ptr;
{

  int f,k,id;

  for (id=0;id<=MAX_ID;id++)
    
    for (f=0;f<FS;f++)
      
      for (k=0;k<=NBBINS;k++)
	{
	  (*bbin_ptr).npackets[k][f][id]=0;
	}
  
  return;
  
}


/*************** find_bbin_intensity***********************/

void find_bbin_intensity(bbin_ptr,lambda_obs_ptr,dlambda_obs_ptr,star, s, iso_flag_ptr)
     
     struct bbin_type *bbin_ptr;
     double *lambda_obs_ptr;
     double *dlambda_obs_ptr;
     struct star_type star[];
     int s;
     int *iso_flag_ptr;
{
  int i,k,f;
  double b,db,dlambda,dlambda2;
  int id;
  
  
   if (MORE==1 || (*iso_flag_ptr)==1 || (RESTART==1 && (*iso_flag_ptr)==0) )
    {
    read_bbins(bbin_ptr);
    }
  

  /* in MJy */
  dlambda=log10(SL2/SL1)/NFBINS;
  for (i =0; i<=NBBINS-1; i++)
    {
      b=((*bbin_ptr).b[i]+(*bbin_ptr).b[i+1])/2.;
      db=(*bbin_ptr).b[i+1]-(*bbin_ptr).b[i];
      for (id=0;id<=MAX_ID;id++)
	for (f=0;f<FS;f++)
	  {
	    /* dlambda2=((pow(10,dlambda/2.)-pow(10,-dlambda/2.))*(*(lambda_obs_ptr+f))); */
	    dlambda2=(*(dlambda_obs_ptr+f));
	    (*bbin_ptr).intensity[i][f][id]=(*bbin_ptr).intensity[i][f][id]+1e17*1e-8*((*(lambda_obs_ptr+f))*(*(lambda_obs_ptr+f))/c_CONST)*(*bbin_ptr).npackets[i][f][id]*(star[s].lum/star[s].npackets)/(2.*M_PI*b*db*4.*M_PI*dlambda2*1e-4);
	  }
    }
  
  (*iso_flag_ptr)=1;
  
  return;
}




/******* make_isogrid *******************************************************/

#ifdef DO_ISOMAP
void make_isogrid(isogrid_ptr)
     
     struct isogrid_type *isogrid_ptr;
     
      
{
  int i,j,f,id,l,m;
    
  (*isogrid_ptr).dx=(XMAX-XMIN)/(1.*NXPIXELS);
  (*isogrid_ptr).dy=(YMAX-YMIN)/(1.*NYPIXELS);
  
/* printf("dx %.3e  dy %.3e\n",dx,dy); */
  

  for (i=0;i<NXPIXELS;i++)
    {
      for  (j=0;j<NYPIXELS;j++)
	
     {
       (*isogrid_ptr).x[i][j]=XMIN+(i+0.5)*(*isogrid_ptr).dx;
       (*isogrid_ptr).y[i][j]=YMAX-(j+0.5)*(*isogrid_ptr).dy;
      
       if (MORE==0)
	 {
	   for (f=0;f<FS;f++)
	     {
	       for (id=0;id<=MAX_ID;id++)
		 {
		   for (l=0;l<NLBINS;l++)
		     for (m=0;m<NMBINS;m++)
		       { (*isogrid_ptr).tau[i][j][l][m]=0.;
			 (*isogrid_ptr).npackets[i][j][f][id][l][m]=0;
			 (*isogrid_ptr).intensity[i][j][f][id][l][m]=0.;
		       }
		 }
	     }
	 }
     }
    }
 
  return;
}


/************** set isogrid.npackets to zero ************************************/
void null_isogrid_npackets(isogrid_ptr)
     
     struct isogrid_type *isogrid_ptr;
     
      
{
  int i,j,f,id,l,m;
  
  for (i=0;i<NXPIXELS;i++)
    for  (j=0;j<NYPIXELS;j++)
      for (f=0;f<FS;f++)
	for (id=0;id<=MAX_ID;id++)
	  for (l=0;l<NLBINS;l++)
	    for (m=0;m<NMBINS;m++)
	      (*isogrid_ptr).npackets[i][j][f][id][l][m]=0;
  
  return;
}

#endif



/*************** function place_into_isogrid ***********************/
/* puts photon  in an observer-dependent grid                      */ 

#ifdef DO_ISOMAP 
void place_into_isogrid(p,i,isogrid_ptr,theta,phi,fbin_obs_ptr, fbin_obs_min_ptr, fbin_obs_max_ptr,l,m,star,s)

     struct photon_packet_type p[];
     int i;
     struct isogrid_type *isogrid_ptr;
     double theta,phi;
     int *fbin_obs_ptr; /* obsrveving lambda bin (if -1 then all) */
     int *fbin_obs_min_ptr;
     int *fbin_obs_max_ptr;
     int l,m;
     struct star_type star[];
     int s;
{
  
  int k,j,id,f;
  double x,y,angle,xx, yy;
  
  int flag1=1;
  int flag2=1;
  double kappa;
  double costheta, cosphi, sintheta, sinphi;
 double theta_new,phi_new;
  
  
  if  (star[s].type==3 && p[i].r==star[s].z)
    
    return;
  
  else
    {
      
      /* calculate the coordinates of the position of the photon before it
	 escapes in the observer's system */
      
      /*
	if (DIM==2) 
	{
	if (p[i] .kz<0) p[i].lz=-p[i].lz; 
	xx=p[i].lx* (p[i].kx/pow(1-p[i].kz*p[i].kz,0.5))+p[i].ly*(p[i].ky/pow(1-p[i].kz*p[i].kz,0.5));
	yy=-p[i]. lx*(p[i].ky/pow(1-p[i].kz*p[i].kz,0.5))+p[i].ly*(p[i].kx/pow(1-p[i].kz*p[i].kz,0.5));
	p[i].lx=xx ;
	p[i].ly=yy; 
	}	
	sinphi=(p[i].ky/pow(1-p[i].kz*p[i].kz,0.5)); 
	cosphi=(p[i].kx/pow(1-p[i].kz*p[i].kz,0.5)); 
	
	sinphi=0.;
	cosphi=1.;
	sintheta=pow(1-p[i].kz*p[i].kz,0.5);
	if (p[i].kz<0) costheta=-p[i].kz; else costheta=p[i].kz;
      */

      /* 
	 theta_new=acos(p[i].kz);
	 
	 cosphi=p[i].kx/sin(theta_new);
	 sinphi=p[i].ky/sin(theta_new);        
	 sintheta=sin(theta_new);  
	 costheta=cos(theta_new); 
      */

      /*
	printf("+kx:% .2f ky:% .2f  kz:% .2f\n",p[i].kx-sin(theta)*cos(phi),p[i].ky-sin(theta)*sin(phi),p[i].kz-cos(theta)); 
	printf("theta1:% .2f phi1:% .2f\n",sin(theta)-sin(theta_new),sin(phi)-sinphi);
	prin`tf("theta-:% .2f phi-:% .2f\n",theta, phi); 
      */

      /* use theta, phi obs cause they define a unique grid. if theta=0 dtheta=10
	 then phi can be 0-->2pi -->everything is rotating around */
      
      /* works only if there is a small photon collection angle */
      x=-sin(phi)*p[i].x+cos(phi)*p[i].y;     
      y=-cos(theta)*cos(phi)*p[i].x-cos(theta)*sin(phi)*p[i].y+sin(theta)*p[i].z;  

       
     

      
 	x=-sin(phi)*p[i].lx+cos(phi)*p[i].ly;     
	y=-cos(theta)*cos(phi)*p[i].lx-cos(theta)*sin(phi)*p[i].ly+sin(theta)*p[i].lz;   
      
     
     

      /* -ant's version ---strange that is uses the direction of the photon and not 
	 the position of the last interaction
	 -the direction of the photon is more or less the direction of the observer 
	 (i am inside the if condition for counting photons
	 -may need modification in 3D problem to get the wanted orientation.
      */
      
     /*  works alright though  */

	kappa=pow(p[i].kx*p[i].kx+p[i].ky*p[i].ky,0.5); 
	
	x=(-p[i].ky*p[i].x+p[i].kx*p[i].y)/kappa;    
	y=(-p[i].kx*p[i].kz*p[i].x-p[i].kz*p[i].ky*p[i].y+kappa*kappa*p[i].z)/kappa;  
      

      flag1=0;
      flag2=0;
      
      
      if (!(x>XMAX || y>YMAX || x<XMIN || y<YMIN))
	
	{
	  
	  k=binary_search_isogrid_x(x,0,NXPIXELS-1,isogrid_ptr,(*isogrid_ptr).dx);
	  j=binary_search_isogrid_y(y,0,NYPIXELS-1,isogrid_ptr,(*isogrid_ptr).dy);
	  
	  for (f=0;f<FS;f++)
	    {
	      
	      if ( (p[i].n>=*(fbin_obs_min_ptr+f) &&  p[i].n<=*(fbin_obs_max_ptr+f)) ||   (*(fbin_obs_ptr+f))<0)   
		{
		  id=p[i].id;
		  
		  (*isogrid_ptr).npackets[k][j][f][id][l][m]=(*isogrid_ptr).npackets[k][j][f][id][l][m]+p[i].weight;	    
		 
		  if (DIM==2 || ISODIM==2)
		    { 
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
			(*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight; /* left-right symmetry */

		      if ((theta==0) || fabs(theta-(M_PI/2.))<1e-5)
			{ 
			  (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			    (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight; /* up-down symmetry */     
			  
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			(*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
			}
		    }	
		  
		  if (DIM==4 || ISODIM==4)
		    { 
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight;
		      
		      if ((theta==0) || fabs(theta-M_PI)<1e-5)
			{ 
			  
			  (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			    (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight;   
		      
			  (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			    (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
			}
		    }
		 
		id=0; /* all packets */
	
		(*isogrid_ptr).npackets[k][j][f][id][l][m]=(*isogrid_ptr).npackets[k][j][f][id][l][m]+p[i].weight;	    
		
		if (DIM==2 || ISODIM==2)
		  { 
		    (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight; /* left-right symmetry */
		    
		    if ((theta==0) || fabs(theta-(M_PI/2.))<1e-5)
		      { 
			(*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			  (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight; /* up-down symmetry */     
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			(*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
		      }
		  }	
		
		if (DIM==4 || ISODIM==4 )
		  { 
		    (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
		      (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight;
		
		    if ((theta==0) || fabs(theta-M_PI)<1e-5)
		      { 
			
			(*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			  (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight;      
			(*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			  (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
		      }
		  }
		
     if (!(p[i].id==3 && (star[s].type==1 || star[s].type==2) && (ISRF_IGNORE_SCATTER==1)))
	
		if (p[i].origin==1) /* just the radiation coming from a specific region in the system */
		  
		  { 
		    id=2;
		    
		    (*isogrid_ptr).npackets[k][j][f][id][l][m]=(*isogrid_ptr).npackets[k][j][f][id][l][m]+p[i].weight;
		    
		    if (DIM==2 || ISODIM==2)
		      { 
			(*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
			  (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight; /* left-right symmetry */
 			
			if ((theta==0) || fabs(theta-(M_PI/2.))<1e-5)
			  { 
			    (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			      (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight; /* up-down symmetry */ 
			    
			    (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			      (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][2][l][m]+p[i].weight; 
			  }
		      }
		    
		    if (DIM==4 || ISODIM==4 )
		      { 
			(*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]=
			  (*isogrid_ptr).npackets[NXPIXELS-k-1][j][f][id][l][m]+p[i].weight;
			
			if ((theta==0) || fabs(theta-M_PI)<1e-5)
			  { 
			    
			    (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]=
			      (*isogrid_ptr).npackets[k][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
			    
			    (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]=
			      (*isogrid_ptr).npackets[NXPIXELS-k-1][NYPIXELS-j-1][f][id][l][m]+p[i].weight; 
			  }
		      }
		  }
		}
	    }
	}
      
      return;
      
    }
}


	   
#endif
/******************** modified binary search **************************/
/********* binary search isogrid y ************************************/
#ifdef DO_ISOMAP
int binary_search_isogrid_x(x, low, high,isogrid_ptr,dx)

     double x;
     int low;
     int high;
     struct isogrid_type *isogrid_ptr;
     double dx;
     
{

  int middle;

  while (low<=high) {

    middle= (low+high)/2;

    if (x>=(*isogrid_ptr).x[middle][0]-dx/2. && x<=(*isogrid_ptr).x[middle][0]+dx/2.) 
      
      return middle;
    
    else if (x<(*isogrid_ptr).x[middle][0]-dx/2.)
      
      high=middle-1;
    
    else low=middle+1;
    
  }
  /* printf("x=%.4e\n",x);
  printf("SEARCH KEY NOT FOUND!\n");
  exit(0);*/

}

/********* binary search isogrid y ******************************/

int binary_search_isogrid_y(y, low, high,isogrid_ptr,dy)

     double y;
     int low;
     int high;
     struct isogrid_type *isogrid_ptr;
     double dy;
{

  int middle;
  
  while (low<=high) 
    {
      middle= (low+high)/2; /* printf("%.3e %d %.3e %.3e\n ",y, middle,(*isogrid_ptr).y[0][49],(*isogrid_ptr).y[0][0]); */
    
      if (y>=(*isogrid_ptr).y[0][middle]-dy/2. && y<=(*isogrid_ptr).y[0][middle]+dy/2.) 

	return middle;

      else if (y<(*isogrid_ptr).y[0][middle]-dy/2.)

	low=middle+1;

      else high=middle-1;

    }
  /*  printf("%.3e %d %.3e %.3e\n ",y, middle,(*isogrid_ptr).y[0][49],(*isogrid_ptr).y[0][0]);
  printf("y=%.4e\n",y);
  printf("SEARCH KEY NOT FOUND2!\n");
  exit(0); */
  
}
#endif



/******* find_isogrid_intensity *******************************************************/


#ifdef DO_ISOMAP
void find_isogrid_intensity(isogrid_ptr,lambda_obs_ptr, background_ptr,dbin,dlambda_obs_ptr,star,s,r_max,iso_flag_ptr,startype)
     
     struct isogrid_type *isogrid_ptr;
     double *lambda_obs_ptr;
     double *background_ptr;	
     struct dbin_type dbin;
     double *dlambda_obs_ptr;    
     struct star_type star[];
     int s;
     double r_max;
     int *iso_flag_ptr;

{
  int i,j,f,id,k,index;int l,m;
  double ds,dlambda ,dlambda2, dr,b,db,r, npixels[NBBINS+1][NFBINS+1][MAX_ID+1],domega,theta,phi;
  struct bbin_type bbin;
 
  if (MORE==1 || (*iso_flag_ptr)==1 || (RESTART==1 &&(*iso_flag_ptr)==0) )
    {
      for (l=0;l<NLBINS;l++)
	for (m=0;m<NMBINS;m++)
	  { 
	    read_isogrid(isogrid_ptr,l,m); 
	  }
    }
  

for (l=0;l<NLBINS;l++)
   for (m=0;m<NMBINS;m++)
     {
       
       domega=dbin.domega[l][m];

#ifdef SPHERICAL_ISOMAP
       domega=4*M_PI;
#endif 
       theta=dbin.theta[l][m];
       phi=dbin.phi[l][m];

       ds=(*isogrid_ptr).dx*(*isogrid_ptr).dy;
       dlambda=log10(SL2/SL1)/NFBINS;
       
       for (i=0;i<NXPIXELS;i++)  
	 
	 for  (j=0;j<NYPIXELS;j++)
	   {
	     for (id=0;id<=MAX_ID;id++)
	       
	       for (f=0;f<FS;f++)
		 {
		   if (DIM==2 || ISODIM==2) 
		     {
		       if ((theta==0) || fabs(theta-(M_PI/2.))<1e-5) 
			 (*isogrid_ptr).npackets[i][j][f][id][l][m]=(*isogrid_ptr).npackets[i][j][f][id][l][m]/4.; 
		       else 
			 (*isogrid_ptr).npackets[i][j][f][id][l][m]=(*isogrid_ptr).npackets[i][j][f][id][l][m]/2.;
		     }  
		   
		   if (DIM==4 || ISODIM==4 ) 
		     {
		       if ((theta==0) || fabs(theta-M_PI)<1e-5) 
			 (*isogrid_ptr).npackets[i][j][f][id][l][m]=(*isogrid_ptr).npackets[i][j][f][id][l][m]/4.; 
		       else 
			 (*isogrid_ptr).npackets[i][j][f][id][l][m]=(*isogrid_ptr).npackets[i][j][f][id][l][m]/2.;
		     }  
		   
		   dlambda2=*(dlambda_obs_ptr+f);
		  
		   (*isogrid_ptr).intensity[i][j][f][id][l][m]=(*isogrid_ptr).intensity[i][j][f][id][l][m]
		     +1e17*1e-8*((*(lambda_obs_ptr+f))*(*(lambda_obs_ptr+f))/c_CONST)*
		     (*isogrid_ptr).npackets[i][j][f][id][l][m]*(star[s].lum/star[s].npackets)/(ds*dlambda2*1e-4*domega);

		   /* if (id==4 && (*isogrid_ptr).npackets[i][j][f][id][l][m]!=0. ) 
		     printf("%.3e %.3e %.3e %.3e %.3e %d \n", (*isogrid_ptr).npackets[i][j][f][id][l][m],
		     (*isogrid_ptr).intensity[i][j][f][id][l][m],ds,dlambda2,domega, f);*/
		     /* I_nu */
		 }
	   }
       
       
       
       /* in case i look from an angle that i expect spherical symmetry:
	  
       i.  devide into bbins
       ii. caclulate how many pixels are in a specific bbin 
       iii.calculate the average intensity of those pixels
       iv. change the intensity of those pixels to the average intensity
       */ 
     
       
       if ((((DIM==2) || (DIM==4)) && ((theta==0) || fabs(theta-M_PI)<1e-5) && EXTRA_SYMMETRY==1) || SPHERICAL_ISOMAP==1)
	 
	 {
	
	 
	   
	   bbin.b[0]=0.;
	   
	   dr=1.*r_max/NBBINS;
	   
	   for(k=1;k<=NBBINS;k++) bbin.b[k]=k*dr; 
	   
	   bbin.b[NBBINS]=r_max;
	   
	   for (id=0;id<=MAX_ID;id++)
	     for (f=0;f<FS;f++) 
	       for (k=0;k<=NBBINS;k++)
		 {
		   bbin.lum[k][f][id]=0.;
		   bbin.npackets[k][f][id]=0;
		   npixels[k][f][id]=0;
		   bbin.intensity[k][f][id]=0.;
		 }
	   
	   for (i=0;i<NXPIXELS;i++)  
	     for  (j=0;j<NYPIXELS;j++)
	       
	       {
		 
		 r=sqrt(pow((*isogrid_ptr).x[i][j],2)+pow((*isogrid_ptr).y[i][j],2));
		 
		 for (id=0;id<=MAX_ID;id++)
		   
		   for (f=0;f<FS;f++)
		     
		for (k=0;k<=NBBINS-1;k++)
		  
		  if (r>=bbin.b[k]-0.5*(*isogrid_ptr).dx && r<bbin.b[k+1]+0.5*(*isogrid_ptr).dx)
		    
		    {
		      bbin.intensity[k][f][id]=bbin.intensity[k][f][id]+(*isogrid_ptr).intensity[i][j][f][id][l][m];
		      npixels[k][f][id]=npixels[k][f][id]+1.;
		      
		    }
	       }
	   
	   for (k=0;k<=NBBINS-1;k++)
	     
	     {
	       /* b=(bbin.b[k]+bbin.b[k+1])/2.;
		  db=bbin.b[k+1]-bbin.b[k]; */
	       
	       for (id=0;id<=MAX_ID;id++)
		 
		 for (f=0;f<FS;f++)
		   {
		     
		     bbin.intensity[k][f][id]=bbin.intensity[k][f][id]/npixels[k][f][id];
		   }
	     }
	   
	   for (i=0;i<NXPIXELS;i++)  
	     
	     for  (j=0;j<NYPIXELS;j++)
	       {
		 r=sqrt(pow((*isogrid_ptr).x[i][j],2)+pow((*isogrid_ptr).y[i][j],2));
		 
		 for (id=0;id<=MAX_ID;id++)
		   
		   for (f=0;f<FS;f++)
		     
		     for (k=0;k<=NBBINS-1;k++)
		       
		       if (r>=bbin.b[k] && r<bbin.b[k+1])
			 
			 (*isogrid_ptr).intensity[i][j][f][id][l][m]=bbin.intensity[k][f][id];
	       }
	 }
      

  /* fix background radiation */
       
       for (i=0;i<NXPIXELS;i++)  
	 for  (j=0;j<NYPIXELS;j++)
	   
	   {
	     for (id=0;id<=1;id++)
	       
	       for (f=0;f<FS;f++)
		 {
		   if (id==3 || id==4)  (*isogrid_ptr).intensity[i][j][f][id][l][m]=0.;
		   else
		     {  
		       if (startype==1 || startype==2 || startype==4) 
			 {

#ifdef CYLINDER_2D 
  if (fabs((*isogrid_ptr).x[i][j])>2*R_MAX_CORE || fabs((*isogrid_ptr).y[i][j])>2*R_MAX_CORE_Z)
#else
 
			   if ( (sqrt(pow((*isogrid_ptr).x[i][j]-(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]-(*isogrid_ptr).dy/2.,2))>=r_max)
				||
				(sqrt(pow((*isogrid_ptr).x[i][j]-(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]+(*isogrid_ptr).dy/2.,2))>=r_max)
				||
				(sqrt(pow((*isogrid_ptr).x[i][j]+(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]-(*isogrid_ptr).dy/2.,2))>=r_max)
				||
				(sqrt(pow((*isogrid_ptr).x[i][j]+(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]+(*isogrid_ptr).dy/2.,2))>=r_max)
				)
			     
#endif
			     (*isogrid_ptr).intensity[i][j][f][id][l][m]=*(background_ptr+f);
			   
			 }
		       
		       if (star[s].type==2) 
			 {
			   
			   if ((sqrt(pow((*isogrid_ptr).x[i][j]-(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]-(*isogrid_ptr).dy/2.,2))>=r_max)
			       || (sqrt(pow((*isogrid_ptr).x[i][j]+(*isogrid_ptr).dx/2.,2)+pow((*isogrid_ptr).y[i][j]+(*isogrid_ptr).dy/2.,2))>=r_max))
			     
			     (*isogrid_ptr).intensity[i][j][f][id][l][m]=lblack((*(lambda_obs_ptr+f))*1e-4,STAR_TEMP)*
			       1e-8*((*(lambda_obs_ptr+f))*(*(lambda_obs_ptr+f))/c_CONST)*1e17;
			 }
		     }
		 }
	   }
       
     }
 

 (*iso_flag_ptr)=1;
 return;

}
#endif

/***** function   *****************************************/

/* cgs units */
void write_sed(fbin_ptr,bisrf_ptr,dbin,r_max,sed_flag_ptr,star_type)

     struct fbin_type *fbin_ptr;
     struct btype *bisrf_ptr;
     struct dbin_type dbin;
     double r_max;
     int *sed_flag_ptr;
     int star_type;
{
  int i,id;
  double tot=0.;
  double lfl,lambda;
  char sed_file[20];char sed_file2[20];
  double dlambda;
  double surface;
  double source[NFBINS+1][NLBINS][NMBINS],background2[NFBINS+1][NLBINS][NMBINS],Lratio,sum1,sum2,back;
  double lambda_array[NFBINS+1],lfl_output[NFBINS+1][NLBINS][NMBINS],lfl_input[NFBINS+1][NLBINS][NMBINS];
  double npackets[NFBINS+1][NLBINS][NMBINS];
  char buff[100];
  double domega;
  int l,m;
  FILE *fp;
   

  for (l=0;l<NLBINS;l++)
    for (m=0;m<NMBINS;m++)
      {	 

	domega=dbin.domega[l][m];

	for (id=0;id<=MAX_ID;id++) 
	  {
#ifdef SPHERICAL
	    surface=4*M_PI*pow(DISTANCE,2);
#endif
	    
#ifndef SPHERICAL
	    surface=domega*pow(DISTANCE,2);
#endif

#ifdef SPHERICAL_SED
	    surface=4*M_PI*pow(DISTANCE,2);
#endif 
   

	    sprintf(sed_file,"sed%d.%d.%d.dat",id,l,m);

	    for (i =0; i<=NFBINS; i++)
	      {
		lambda_array[i]=0;
		lfl_output[i][l][m]=0;
		lfl_input[i][l][m]=0;
		npackets[i][l][m]=0;
		source[i][l][m]=0.;
		background2[i][l][m]=0.;
	      }
  
	    
	    if (MORE==1 || (*sed_flag_ptr)==1 || (RESTART==1 && (*sed_flag_ptr)==0) )
	      {
		
		if (( fp = fopen(sed_file, "r")) != NULL)
		  {
		    for (i = 0; (i <=NFBINS) & getdata(fp, buff); i++)
		      
		      sscanf(buff,"%lf %lf %lf ",
			     &lambda_array[i],
			     &lfl_output[i][l][m],
			     &npackets[i][l][m]);
		    
		    fclose(fp);
		  }
		else
		  {

		    fprintf(stderr, "\nWARNING: FAILED TO OPEN %s FILE. \n\n",sed_file); exit(0);
		  }
	      }

	    /* Open  file */
	    if (( fp = fopen(sed_file, "w")) == NULL)
	      {fprintf(stderr, "FAILED TO OPEN %s FILE. \n",sed_file);return;}
	    
	    dlambda=log10(SL2/SL1)/NFBINS;

	    if (star_type==1 || star_type==2 || star_type==4) /*  ISRF */ 
	      {
		   
		fprintf(fp,"@lambda (microns), lum*lambda (erg/s cm2), npackets  \n");
		
		for (i =0; i<=NFBINS; i++)
		  {
		    lambda=(*fbin_ptr).lambda[i];
		    lfl=(*fbin_ptr).lum[i][id][l][m]/((pow(10,dlambda/2.)-pow(10,-dlambda/2.))*lambda);
		    
		    fprintf(fp,"%7.3f %.4e %.4e \n",
			    lambda,lfl_output[i][l][m]+(lfl*lambda/surface),
			    (*fbin_ptr).npackets[i][id][l][m]);
		    tot=tot+(*fbin_ptr).lum[i][id][l][m];
		    
		    if (id==2) 
		      {
			source[i][l][m]=lfl_output[i][l][m]+(lfl*lambda/surface); /* source flux */
			background2[i][l][m]= M_PI*isrf(lambda,bisrf_ptr)*1.e4*c_CONST/lambda*pow(r_max/DISTANCE,2);
		      }
		  }
	      }
	    
	    else
	      
	      {
		/* STAR radiation field */
		  
		fprintf(fp,"@lambda, lum*lambda/F, npackets\n");
		
		for (i =0; i<=NFBINS; i++)   
		  { 
		    
		    lambda=(*fbin_ptr).lambda[i];
		    lfl=(*fbin_ptr).lum[i][id][l][m]/((pow(10,dlambda/2.)-pow(10,-dlambda/2.))*lambda);
		    
		    
		    fprintf(fp,"%7.3f %.4e %.4e \n",
			    (*fbin_ptr).lambda[i],
			    lfl_output[i][l][m]+(lfl*(*fbin_ptr).lambda[i]/surface), (*fbin_ptr).npackets[i][id][l][m]);
		    
		    
		    if (id==2) 
		      {
			source[i][l][m]=lfl_output[i][l][m]+(lfl*lambda/surface); /* source flux */
			background2[i][l][m]= M_PI*isrf(lambda,bisrf_ptr)*1.e4*c_CONST/lambda*pow(r_max/DISTANCE,2);
		      }
		    		    
		  }
		  
	      }
	  
	    
	    /* SPEC OUTPUT 

	       printf("-- l:%d m:%d id:%d -- [solid angle (of 4pi):%.4e %s]\n",l,m,id,100*domega/(4*M_PI),"%");
	       printf("Lum of the photons contributing to SED: %.4e\n", tot);

	    */
	    fclose(fp);
	    
	    /* calculate Lratio and output src file with just the source */
	    if (id==2)
	      {
		sum1=0.;
		sum2=0.;
		back=0.;
		
		for (i =0; i<=NFBINS; i++)
		  {
		    if ((*fbin_ptr).lambda[i]>=350) sum1=fabs(source[i][l][m])*(pow(10,dlambda/2.)-pow(10,-dlambda/2.))*surface+sum1;
		    sum2=fabs(source[i][l][m])*(pow(10,dlambda/2.)-pow(10,-dlambda/2.))*surface+sum2;
		    back=back+background2[i][l][m]*(pow(10,dlambda/2.)-pow(10,-dlambda/2.))*surface;
		  }
		Lratio=100.*sum1/sum2;
		
		sprintf(sed_file2,"src.sed%d.dat",id);
		
		/* Open  file */
		if (( fp = fopen(sed_file2, "w")) == NULL)
		  {fprintf(stderr, "FAILED TO OPEN sed_file2 FILE. \n");return;}
		
		fprintf(fp,"@lambda (microns), lum*lambda (erg/s cm2)(source), npackets, lambda F_lambda (background)(erg/s cm2) \n");
		
		for (i =0; i<=NFBINS; i++)  fprintf(fp,"%7.3f %.4e %6d %.4e \n",(*fbin_ptr).lambda[i], source[i][l][m], 0, background2[i][l][m]);
		
		fprintf(fp,"@ Lsmm/Lbol=%.5f\n",Lratio);
		
		fclose(fp);		/* SPEC OUTPUT	
		printf("Lsmm/Lbol=%.5f, Source Lum=%.3e (Lsun), Background Lum=%.3e (Lsun) \n",
		       Lratio,sum2/(3.862*1e33),back/(3.862*1e33));
		*/
		       }
	    
	  }
      
	
	/*  printf("---------------------------------------------------------------------------\n"); */
}

  (*sed_flag_ptr)=1;

  return;
}

#ifdef POST_RT

/***** function write_sed_prt  *****************************************/

/* post radiative transfer SED */

/* cgs units */
void write_sed_prt(fbin_ptr,bisrf_ptr,dbin,r_max,c,starcell,cextra,star,st,ncells)
  
     struct fbin_type *fbin_ptr;
     struct btype *bisrf_ptr;
     struct dbin_type dbin;
     double r_max;
     struct cell_type c[]; 
     struct cell_type starcell[];
     struct cell_type cextra[]; 
     struct star_type star[];
     int st;
int ncells;
{
  int i,cell_i,j,k,l,m;
  long cex_i,s;
  double tot=0.;
  double lambda;
  char sed_file[10];
  double mass1=0.,mass2=0.;
  double dlambda,domega;
  double surface;
  double source[NFBINS+1][NLBINS][NMBINS],background2[NFBINS+1][NLBINS][NMBINS];
  double lfl_output[NFBINS+1][NLBINS][NMBINS],lfl_input[NFBINS+1][NLBINS][NMBINS];
  int npackets[NFBINS+1][NLBINS][NMBINS];
  int id;
  FILE *fp;
  
  for (l=0;l<NLBINS;l++)
    for (m=0;m<NMBINS;m++)
      {	 
	
	domega=dbin.domega[l][m];
	id=0;

	
	surface=domega*pow(DISTANCE,2); 
	
	sprintf(sed_file,"sed%d.%.0f.%.0f.prt.dat",id,dbin.theta[l][m]*180/M_PI,dbin.phi[l][m]*180/M_PI);
	
	for (i =0; i<=NFBINS; i++)
	  {
	    
	    lfl_output[i][l][m]=0;
	    lfl_input[i][l][m]=0;
	    npackets[i][l][m]=0;
	    source[i][l][m]=0.;
	    background2[i][l][m]=0.;
	  }
	
#ifdef DUST_MIX
	printf(" HEY STOP IT! DUST MIX DOES NOT WORK ... FIX IT\n");
#endif
	
	/* Open  file */
	if (( fp = fopen(sed_file, "w")) == NULL)
	  {fprintf(stderr, "FAILED TO OPEN sed_file FILE. \n");return;}
	
	dlambda=log10(SL2/SL1)/NFBINS;
	
	fprintf(fp,"@lambda, lum*lambda/F, npackets, lambda F_lambda/F (input))\n");
  
	for (i =0; i<=NFBINS; i++)   
	  {       
	    lambda=(*fbin_ptr).lambda[i];
	    /* lfl=(*fbin_ptr).lum[i][id][l][m]/((pow(10,dlambda/2.)-pow(10,-dlambda/2.))*lambda); */
	    
	    if (DO_SPH_CELLS==1)	  /* normal sph cells */ 
	      {  
		for (cell_i=0;cell_i<ncells-star[st].rcells*star[st].lcells;cell_i++)
		  { 
		    if (c[cell_i].temp>0 && 
			(c[cell_i].x-star[st].x)*(c[cell_i].x-star[st].x)+
			(c[cell_i].y-star[st].y)*(c[cell_i].y-star[st].y)+
			(c[cell_i].z-star[st].z)*(c[cell_i].z-star[st].z)<=ORIGIN_SIZE*ORIGIN_SIZE)
		      {
			for(s=0;s<8;s++)
			  {
			    lfl_output[i][l][m]= lfl_output[i][l][m]+(*fbin_ptr).lambda[i]*1e-4*
			      (*fbin_ptr).abs[i]*(c[cell_i].mass/8.)*
			      lblack((*fbin_ptr).lambda[i]*1e-4,c[cell_i].temp)*
			      exp(-(*fbin_ptr).opacity[i]*c[cell_i].prt_tau[l][m][s])/(DISTANCE*DISTANCE); 
			  }
		      }
		  }
	      }	  
	    
	    /*	    printf(" rmax_2/(DISTANCE*DISTANCE) ????\n");*/
	    
	    if (DO_STAR_CELLS==1) /* star.grid cells  */
	      {
		mass2=0.;
		for (k=ncells-star[st].rcells*star[st].lcells;k<ncells;k++)             
		  {      
		    
		    if (c[k].temp>0) mass1=mass1+c[k].mass;
		    
		    for (j=0;j<STAR_MCELLS;j++)
		      {
			cex_i=j+(k-(ncells-star[st].rcells*star[st].lcells))*star[st].mcells;
			if (cextra[cex_i].temp>0)
			  { 
			    mass2=mass2+cextra[cex_i].mass;
			    
			    lfl_output[i][l][m]= lfl_output[i][l][m]+(*fbin_ptr).lambda[i]*1e-4*
			      (*fbin_ptr).abs[i]*cextra[cex_i].mass*
			      lblack((*fbin_ptr).lambda[i]*1e-4,cextra[cex_i].temp)*
			      exp(-(*fbin_ptr).opacity[i]*(cextra[cex_i].prt_tau[l][m][0]))/(DISTANCE*DISTANCE); 
			  }	
		      }
		  } 
	      }
	    
	    /* printf("mass2=%.4e\n",mass2/SOLAR_MASS); */
	    
	    if (DO_STAR==1) /* star flux */
	      {
		lfl_output[i][l][m]=lfl_output[i][l][m]+
		  (*fbin_ptr).lambda[i]*1e-4*
		  lblack((*fbin_ptr).lambda[i]*1e-4,star[st].temp)*
		  exp(-(*fbin_ptr).opacity[i]*starcell[0].prt_tau[l][m][0])*
		  (4*M_PI*star[st].radius*star[st].radius)/(4*M_PI*DISTANCE*DISTANCE)*M_PI;
		
		
	      }
	    
	    fprintf(fp,"%7.3f %.4e %.4e %.4e \n",
		    (*fbin_ptr).lambda[i],
		    lfl_output[i][l][m], 0.,
		    star[st].dilution*M_PI*((domega*star[st].radius*star[st].radius)/surface)
		    *((*fbin_ptr).lambda[i]*1e-4)*lblack((*fbin_ptr).lambda[i]*1e-4,1.*star[st].temp));
	    
	    tot=tot+(*fbin_ptr).lum[i][id][l][m];	      
	    
	  }
	
      	fclose(fp);
      }
 
  return;
  
} 

#endif

#ifdef DO_ISOMAP


/********** write isogrid for IDL **************************************************/

void isomap_idl(isogrid_ptr)

 struct isogrid_type *isogrid_ptr;

{

  int i,j,l,k,m,id;                       
  char      filename2[20];

  FILE *fp;


  
  #ifdef SPHERICAL
  
  /**** write isophotal map (ready for IDL) ****/
	l=0;
	m=0;
	
	for (id=0;id<=MAX_ID;id++)
	  {
	    for (k=0;k<FS;k++)
	      {
		sprintf(filename2,"%s.%d.%d.%d.%d.dat","iso",k,id,l,m);	       
		if (( fp = fopen(filename2, "w")) == NULL)
		  {fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename2);exit(0);}	       
	       for (i=NYPIXELS-1 ; i>=0; i--)
		 {
		   for (j =0; j<NXPIXELS; j++) 
		     {		       		  
		       fprintf(fp,"% .3e",(*isogrid_ptr).intensity[j][i][k][id][l][m]);
		     }
		   fprintf(fp,"\n");
		 }	       
	       fclose(fp);
	     }	
  }
  #else
  for (l=0;l<NLBINS;l++)
    
    for (m=0;m<NMBINS;m++)
      { 
	
	/**** write isophotal map (ready for IDL) ****/
	
	for (id=0;id<=MAX_ID;id++)
	  {
	    for (k=0;k<FS;k++)
	      {
		sprintf(filename2,"%s.%d.%d.%d.%d.dat","iso",k,id,l,m);	       
		if (( fp = fopen(filename2, "w")) == NULL)
		  {fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename2);exit(0);}	       
	       for (i=NYPIXELS-1 ; i>=0; i--)
		 {
		   for (j =0; j<NXPIXELS; j++) 
		     {		       		  
		       fprintf(fp,"% .3e",(*isogrid_ptr).intensity[j][i][k][id][l][m]);
		     }
		   fprintf(fp,"\n");
		 }	       
	       fclose(fp);
	     }	
	 }       
       
      }
#endif




  return;
}

/********** write statistical noise files  for IDL ********************************/

void isomap_statistics_idl(isogrid_ptr,s)
     
     struct isogrid_type *isogrid_ptr;
     int s;

{/*** statistical noise files (not good for spherical case!!!) ***/
   int i,j,k,l,m,id;                       
  char      filename2[20];

  FILE *fp;
  
 
  
  for (l=0;l<NLBINS;l++)
    
    for (m=0;m<NMBINS;m++)
      { 
	
       /**** write isophotal map (ready for IDL) ****/
	
       for (id=0;id<=MAX_ID;id++)
	 {
	   for (k=0;k<FS;k++)
	     {
	       sprintf(filename2,".%d%s.%d.%d.%d.%d.dat",s,"stat.star",k,id,l,m);
	       
	       if (( fp = fopen(filename2, "w")) == NULL)
		 {
		   fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename2);
		   exit(0);
		 }
	       
	       for (i=NYPIXELS-1 ; i>=0; i--)
		 {
		   for (j =0; j<NXPIXELS; j++) 
		     {		       
		       if ((*isogrid_ptr).npackets[j][i][k][id][l][m]>0.)
			 fprintf(fp,"%.3f\n",100./sqrt((*isogrid_ptr).npackets[j][i][k][id][l][m]));
		       else
			 {
			   if ((*isogrid_ptr).intensity[j][i][k][id][l][m]==(*isogrid_ptr).intensity[0][0][k][id][l][m])
			     fprintf(fp,"0.\n");
			   else
			     fprintf(fp,"100.\n");
			 }
		       
		       
		     }
		   fprintf(fp,"\n");
		 }	       
	       fclose(fp);
	     }
	   
	 }          
      }

  return;
}
  
#endif


/********** calculate domega of the photons observed **********************/

void calculate_domega(dbin_ptr)
     struct dbin_type *dbin_ptr;
{
	  	
  int l,m;
  double theta_observer;
  
  if (DIM==3 && ISODIM==3 )
    {
      
      printf("\nGeneral case... 3\n");	 	     
      
      /* solid angle for which photons are counted*/
      
      for (l=0;l<NLBINS;l++)
	for (m=0;m<NMBINS;m++)
	  {
	    (*dbin_ptr).domega[l][m]=2*M_PI*(1-cos(DTHETA)); 
	    printf("omega: %.4e\n",(*dbin_ptr).domega[l][m]);
	    /*this is only when there in only 1 obs_theta */
	  }
    }
  
  if (DIM==2 || ISODIM==2 )
   {	      
     printf("\nGeneral case: phi symmetry, theta -theta symmetry\n");	      
     
     for (l=0;l<NLBINS;l++)
       for (m=0;m<NMBINS;m++)
	 {	       
	   theta_observer=(*dbin_ptr).theta[l][m];    
	   /* problem when i have only one angle =0 and MPI tolerance */
	   
	   if (theta_observer-DTHETA<0) 
	    (*dbin_ptr).domega[l][m]=2*2*M_PI*(1-cos(theta_observer+DTHETA));
	   
	   else if (theta_observer+DTHETA>M_PI) 
	     
	     (*dbin_ptr).domega[l][m]=2*2*M_PI*(cos(theta_observer-DTHETA)+1);
	   
	   else if (fabs(theta_observer-M_PI/2.)<=DTHETA)
	     {
	       if (theta_observer>=M_PI/2.) 
		 
		 (*dbin_ptr).domega[l][m]=2*M_PI*((-2)*cos(theta_observer+DTHETA));
	       
	       else 
		 
		 (*dbin_ptr).domega[l][m]=2*M_PI*(2*cos(theta_observer-DTHETA));
	     }

	   else 
	     
	     (*dbin_ptr).domega[l][m]=2*2*M_PI*
	       (cos(theta_observer-DTHETA)-cos(theta_observer+DTHETA));

printf("omega: %.4e\n",(*dbin_ptr).domega[l][m]);

	 }
   }

  if (DIM==4 || ISODIM==4)
    {
      printf("\nGeneral case: phi symmetry\n");
      
      for (l=0;l<NLBINS;l++)
       for (m=0;m<NMBINS;m++)
	 {	      
	   theta_observer=(*dbin_ptr).theta[l][m];
	   
	   if (theta_observer-DTHETA<0) 
	     
	    (*dbin_ptr).domega[l][m]=2*M_PI*(1-cos(theta_observer+DTHETA));
	   
	   else 
	     {
	       if (theta_observer+DTHETA>M_PI) 
		 
		 (*dbin_ptr).domega[l][m]=2*M_PI*(cos(theta_observer-DTHETA)+1);
	       
	       else 
		 
		 (*dbin_ptr).domega[l][m]=2*M_PI*
		   (cos(theta_observer-DTHETA)-cos(theta_observer+DTHETA));
	     }
	 }
    }
  
  return;
}

/********** find out in which fbins are the given observation wavelenths ****************/
 
void find_obs_fbins(fbin_obs_ptr,fbin, dlambda_obs_ptr,fbin_obs_min_ptr, fbin_obs_max_ptr)
  
     int *fbin_obs_ptr;
     struct fbin_type fbin;
     double *dlambda_obs_ptr;  
     int *fbin_obs_min_ptr;
     int *fbin_obs_max_ptr;
{
  int f,k,i;
  double dlambda;
  double lambda_obs[FS]; // wavelengths to construct isophotal maps 
  double dl_obs[FS]; // wavelengths to construct isophotal maps 

 for (i=0;i<FS;i++)
     {
         switch (i)
           {
           case 0:lambda_obs[0]=F0;dl_obs[0]=dF0/2.;break;
           case 1:lambda_obs[1]=F1;dl_obs[1]=dF1/2.;break;
           case 2:lambda_obs[2]=F2;dl_obs[2]=dF2/2.;break;
           case 3:lambda_obs[3]=F3;dl_obs[3]=dF3/2.;break;
           case 4:lambda_obs[4]=F4;dl_obs[4]=dF4/2.;break;
           case 5:lambda_obs[5]=F5;dl_obs[5]=dF5/2.;break;
           case 6:lambda_obs[6]=F6;dl_obs[6]=dF6/2.;break;
           case 7:lambda_obs[7]=F7;dl_obs[7]=dF7/2.;break;
           case 8:lambda_obs[8]=F8;dl_obs[8]=dF8/2.;break;
           case 9:lambda_obs[9]=F9;dl_obs[9]=dF9/2.;break;
           }

       }
 
  (*(fbin_obs_ptr+0))=-1;
  
  dlambda=log10(SL2/SL1)/NFBINS;
 

printf("\n---------- Wavelength BINS  ---------------------------\n");

    printf("::Details of the wavelegth bins (l1,l2)\n");
 
  for (f=0;f<FS;f++)
    
    for (k=0;k<=NFBINS;k++)       
      
      if (lambda_obs[f]>(pow(10,-dlambda/2.)*fbin.lambda[k]) && 
	  (lambda_obs[f]<=pow(10,dlambda/2.)*fbin.lambda[k]) )
	
	{
	  *(fbin_obs_ptr+f)=k;
	  

	  (*(dlambda_obs_ptr+f))=fbin.lambda2[k]-fbin.lambda1[k];

	  printf("\nf:%d FBIN=%d...lambda1=%.3f lambda2=%.3f dlambda=%.3f\n",
		 f,*(fbin_obs_ptr+f),fbin.lambda1[k], fbin.lambda2[k], *(dlambda_obs_ptr+f));
	  break;
	}
  

  for (f=0;f<FS;f++)

    {
 *(fbin_obs_min_ptr+f)=10000000;
 *(fbin_obs_max_ptr+f)=0;

    for (k=0;k<=NFBINS;k++)       
      {
      if (lambda_obs[f]+dl_obs[f]>(pow(10,-dlambda/2.)*fbin.lambda[k]) && 
	  (lambda_obs[f]-dl_obs[f]<=pow(10,dlambda/2.)*fbin.lambda[k]) )
	
	{
	  if (k<=*(fbin_obs_min_ptr+f))  *(fbin_obs_min_ptr+f)=k;
	  if (k>=*(fbin_obs_max_ptr+f))  *(fbin_obs_max_ptr+f)=k;
	}
      }
     
      (*(dlambda_obs_ptr+f))=fbin.lambda2[*(fbin_obs_max_ptr+f)]-fbin.lambda1[*(fbin_obs_min_ptr+f)];

      printf("\nf:%d  fbin_min:%d  fbin_max:%d  dlambda:%.3f   (%.3f-%.3f micron)\n" ,f, *(fbin_obs_min_ptr+f),*(fbin_obs_max_ptr+f),(*(dlambda_obs_ptr+f)),fbin.lambda1[*(fbin_obs_min_ptr+f)], fbin.lambda2[*(fbin_obs_max_ptr+f)]  );
     }
  return;
}

/***********************************************************************************/

void find_fbin_isogrid_3D(p,i,dbin,isogrid_ptr,fbin_ptr,fbin_obs_ptr, fbin_obs_min_ptr, fbin_obs_max_ptr,total_fbin_ptr,star,s)
  
     struct photon_packet_type p[];
int i;
struct dbin_type dbin;
struct isogrid_type *isogrid_ptr;
struct fbin_type *fbin_ptr;
int *fbin_obs_ptr; /* observing lambda bin (if -1 then all) */
int *fbin_obs_max_ptr; 
int *fbin_obs_min_ptr; 
double *total_fbin_ptr;
struct star_type star[];
int s;

{
  int l,m;
  double phi_observer;
  double theta_observer;
  double angle;
  double cosorigintheta;

  for (l=0;l<NLBINS;l++)
    for (m=0;m<NMBINS;m++)
      {
	
	/* printf("theta:%.2f phi:%.2f\n",dbin.theta[l][m]*180/M_PI,dbin.phi[l][m]*180/M_PI); */
	
	phi_observer=dbin.phi[l][m];
	theta_observer=dbin.theta[l][m];
	angle=acos(cos(phi_observer)*sin(theta_observer)*p[i].kx+
		   sin(phi_observer)*sin(theta_observer)*p[i].ky+
		   cos(theta_observer)*p[i].kz); 
	
	/* printf("ANGLE:%.2f\n",angle*180/M_PI); */   
#ifdef SPHERICAL_SED
#ifdef ORIGIN  /* takes a circural area on the projected image and calculates EVERYTHING that comes inside there */
               /* there is a problem!!!!!!!!!!!!!!!! */
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    p[i].origin=0;
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif	    

	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);
#endif

#ifdef  SPHERICAL_ISOMAP
#ifdef  DO_ISOMAP
  	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,l,m,star,s);
#endif
#endif	
	if (angle<=DTHETA) 
	  {	   

#ifndef SPHERICAL_SED
#ifdef ORIGIN  /* takes a circural area on the projected image and calculates EVERYTHING that comes inside there */
               /* there is a problem!!!!!!!!!!!!!!!! */
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    p[i].origin=0;
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif	    

	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);
#endif

#ifndef  SPHERICAL_ISOMAP
#ifdef  DO_ISOMAP
  	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,l,m,star,s);
#endif
#endif
	  }
      }
  
  return;
}

/***********************************************************************************/

/* some symmetry: phi symmetry and theta, -theta symmetry */

void find_fbin_isogrid_2D_I(p,i, dbin,isogrid_ptr,fbin_ptr,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,total_fbin_ptr,star,s)
  
     struct photon_packet_type p[];
int i;
struct dbin_type dbin;
struct isogrid_type *isogrid_ptr;
struct fbin_type *fbin_ptr;
int *fbin_obs_ptr; /* observing lambda bin (if -1 then all) */
int *fbin_obs_min_ptr;
int *fbin_obs_max_ptr;
double *total_fbin_ptr;
struct star_type star[];
int s;
{	 
  
  int l,m;
  
  double phi_observer;
  double theta_observer;
  double cosorigintheta;
 
  
  for (l=0;l<NLBINS;l++)
    for (m=0;m<NMBINS;m++)
      {
	phi_observer=dbin.phi[l][m];
	theta_observer=dbin.theta[l][m];

#ifdef SPHERICAL_SED
#ifdef ORIGIN
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    p[i].origin=0;    
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif	    
	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);  
#endif

#ifdef  SPHERICAL_ISOMAP
#ifdef  DO_ISOMAP
  	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,l,m,star,s);
#endif
#endif	

	
	if( (acos(p[i].kz)<=theta_observer+DTHETA && acos(p[i].kz)>=theta_observer-DTHETA)  ||
	    (M_PI-acos(p[i].kz)<=theta_observer+DTHETA && M_PI-acos(p[i].kz)>=theta_observer-DTHETA  ) )       
	  
	  {	  
#ifndef SPHERICAL_SED
#ifdef ORIGIN
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    p[i].origin=0;    
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif	    
	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);       
#endif

#ifndef SPHERICAl_ISOMAP
#ifdef DO_ISOMAP
	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,l,m,star,s);
#endif
#endif	    
	  }
      }
  
  return;
}
/***********************************************************************************/

/* just phi symmetry  */

  void find_fbin_isogrid_2D_II(p,i,dbin,isogrid_ptr,fbin_ptr,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,total_fbin_ptr,star,s)
  
     struct photon_packet_type p[];
int i;
struct dbin_type dbin;
struct isogrid_type *isogrid_ptr;
struct fbin_type *fbin_ptr;
int *fbin_obs_ptr; /* observing lambda bin (if -1 then all) */
int *fbin_obs_min_ptr;
int *fbin_obs_max_ptr;
double *total_fbin_ptr;
struct star_type star[];
int s;
{
  int l,m;
  double phi_observer;
  double theta_observer;
  double cosorigintheta;

  for (l=0;l<NLBINS;l++)
    for (m=0;m<NMBINS;m++)
      {
	phi_observer=dbin.phi[l][m];
	theta_observer=dbin.theta[l][m];

#ifdef SPHERICAL_SED
#ifdef ORIGIN 
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    
	    p[i].origin=0;
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif
	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);      
#endif

#ifdef  SPHERICAL_ISOMAP
#ifdef  DO_ISOMAP
  	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr, fbin_obs_min_ptr,fbin_obs_max_ptr,l,m,star,s);
#endif
#endif	

	if( (acos(p[i].kz)<=theta_observer+DTHETA && acos(p[i].kz)>=theta_observer-DTHETA))
	  {	 

#ifndef SPHERICAL_SED
#ifdef ORIGIN 
	    cosorigintheta=(p[i].lx*sin(dbin.theta[l][m])*cos(dbin.phi[l][m])+
			    p[i].ly*sin(dbin.theta[l][m])*sin(dbin.phi[l][m])+
			    p[i].lz*cos(dbin.theta[l][m]))/p[i].lr;
	    
	    p[i].origin=0;
	    if (fabs(sqrt(1-pow(cosorigintheta,2)))*p[i].lr<ORIGIN) p[i].origin=1;
#endif
	    place_into_fbins(p,i,fbin_ptr,total_fbin_ptr,l,m,star,s);  
#endif
    
#ifndef SPHERICAL_ISOMAP
#ifdef DO_ISOMAP
	    place_into_isogrid(p,i,isogrid_ptr,theta_observer,phi_observer,fbin_obs_ptr,fbin_obs_min_ptr, fbin_obs_max_ptr,l,m,star,s);
#endif
#endif	    
	  }
      }
  return;
}






#ifdef HELL
void find_isogrid_intensity(isogrid_ptr,lambda_obs_ptr, background_ptr,dbin,dlambda_obs_ptr,star,s,r_max)
    
     struct isogrid_type *isogrid_ptr;
     double *lambda_obs_ptr;
     double *background_ptr;	
     struct dbin_type dbin;
     double *dlambda_obs_ptr;    
     struct star_type star[];
     int s;
     double r_max;

{
  int i,j,f,id,k,index;int l,m;double theta,phi;
  double ds,dlambda ,dlambda2, dr,b,db,r, npixels[NBBINS+1][NFBINS+1][MAX_ID+1],domega;
  struct bbin_type bbin;
 

       if ((((DIM==2) || (DIM==4)) && ((theta==0) || fabs(theta-M_PI)<1e-5) && EXTRA_SYMMETRY==1)  || SPHERICAL_ISOMAP==1)
	 
	 {
	   
	  
	   bbin.b[0]=0.;
	   
	   dr=1.*r_max/NBBINS;
	   
	   for(k=1;k<=NBBINS;k++) bbin.b[k]=k*dr; 
	   
	   bbin.b[NBBINS]=r_max;
	   
	   for (id=0;id<=MAX_ID;id++)
	     for (f=0;f<FS;f++) 
	       for (k=0;k<=NBBINS;k++)
		 {
		   bbin.lum[k][f][id]=0.;
		   bbin.npackets[k][f][id]=0;
		   npixels[k][f][id]=0;
		   bbin.intensity[k][f][id]=0.;
		 }
	   
	   for (i=0;i<NXPIXELS;i++)  
	     for  (j=0;j<NYPIXELS;j++)
	       
	       {
		 
		 r=sqrt(pow((*isogrid_ptr).x[i][j],2)+pow((*isogrid_ptr).y[i][j],2));
		 
		 for (id=0;id<=MAX_ID;id++)
		   
		   for (f=0;f<FS;f++)
		     
		for (k=0;k<=NBBINS-1;k++)
		  
		  if (r>=bbin.b[k]-0.5*(*isogrid_ptr).dx && r<bbin.b[k+1]+0.5*(*isogrid_ptr).dx)
		    
		    {
		      bbin.intensity[k][f][id]=bbin.intensity[k][f][id]+(*isogrid_ptr).intensity[i][j][f][id][l][m];
		      npixels[k][f][id]=npixels[k][f][id]+1.;
		      
		    }
	       }
	   
	   for (k=0;k<=NBBINS-1;k++)
	     
	     {
	       /* b=(bbin.b[k]+bbin.b[k+1])/2.;
		  db=bbin.b[k+1]-bbin.b[k]; */
	       
	       for (id=0;id<=MAX_ID;id++)
		 
		 for (f=0;f<FS;f++)
		   {
		     
		     bbin.intensity[k][f][id]=bbin.intensity[k][f][id]/npixels[k][f][id];
		   }
	     }
	   
	   for (i=0;i<NXPIXELS;i++)  
	     
	     for  (j=0;j<NYPIXELS;j++)
	       {
		 r=sqrt(pow((*isogrid_ptr).x[i][j],2)+pow((*isogrid_ptr).y[i][j],2));
		 
		 for (id=0;id<=MAX_ID;id++)
		   
		   for (f=0;f<FS;f++)
		     
		     for (k=0;k<=NBBINS-1;k++)
		       
		       if (r>=bbin.b[k] && r<bbin.b[k+1])
			 
			 (*isogrid_ptr).intensity[i][j][f][id][l][m]=bbin.intensity[k][f][id];
	       }
	 }
      
#endif
