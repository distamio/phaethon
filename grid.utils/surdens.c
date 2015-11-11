/* 

FOR  grids: grid_dim2_c and grid_dim4_c

output files: isodens.dat 
              isotemp.dat 

to be used from IDL

it has to be compild locally

compile: gcc -g -o rdens  surdens.c dens.c /home/athena/spxstd/npha/main/inout.c -lm


*IMPORTANT: there are 2 places to change NXPIXELS in this file if i want
            better pictures. 

*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "/home/athena/spxstd/npha/main/constants.h"
#include "/home/athena/spxstd/npha/main/structs.h"
#include "/home/athena/spxstd/npha/main/inout.h"
#include "/home/athena/spxstd/npha/main/phaethon.functions.h"
#include "/home/athena/spxstd/npha/main/photons.h"
#include "/home/athena/spxstd/npha/main/spec.functions.h"
#include "/home/athena/spxstd/npha/main/functions.h"
#include "/home/athena/spxstd/npha/grid/grid.functions.h"
#include "../main/params.h"
#include "../main/phaethon.params.h"
#include "params.grid.h"
#include "../opacity/params.opa.h"

#define XMIN (-0.8*PC)
#define XMAX  (-XMIN) 
#define YMIN   (XMIN) 
#define YMAX   (XMAX)
#define NXPIXELS 600
#define NYPIXELS 600
#define dx ((XMAX-XMIN)/(1.*NXPIXELS))
#define dy dx
struct isogrid_type2
       {
	 double x[NXPIXELS][NYPIXELS];         
	 double y[NXPIXELS][NYPIXELS];      
         double r[NXPIXELS][NYPIXELS];
         double npackets[NXPIXELS][NYPIXELS][10][5][NLBINS][NMBINS];
         double intensity[NXPIXELS][NYPIXELS][10][5][NLBINS][NMBINS];
       };

void make_isogrid2(struct isogrid_type2 *);



double lambda_obs[FS+1]={F0,F1,F2,F3,F4,F5,F6,F7};

struct btype bisrf;
struct photon_packet_type phtype;
struct cell_type ctype;
struct opacity_type opacity;
struct mean_op_type mean_op;
struct reemit_type reemit;
struct fbin_type fbin;
struct dbin_type dbin;
struct isogrid_type2 isogrid;
long   pdum=-123;      /* initial values for the random number generator           */
long   idum=-1;        /* initial value for random number generator                */
double mfp_free=0;     /* distance the photons travel in a step in  empty space    */
int    ncells=0;       /* number of cells in the grid -1 to run from 0 to NCELLS-1 */
double npackets=0;     /* total.. normal +multiple                                 */


int    interactions=0; /* total number of interactions                             */
double time0=0;        /* starting cpu clock time */
double time1=0;        /* ending cpu clock time   */
double intercount[21];/* percentage of photons that interact 0,1,2...,20 times    */


struct fbin_type fbin;
struct dbin_type dbin;
struct bbin_type bbin;


int       ndirect=0;          /* number of direct photons escaping */
double    total_fbin=0.;      /* number of photons placed in frequency bins     */
double    total_dbin;     /* number of photons placed in direction bins     */
double    total_bbin[10];     /* number of photons placed in b bins             */

double    lbin_obs=0;  
double    mbin_obs=0;  
int       fbin_obs[10];

char      filename[20],filename1[20],filename2[20],filename3[20];
double    theta_observer;
double    phi_observer;

double    background[10];

int main()

{
double photon_i=0,dlambda,r,theta,c1,c2;         /* photon counter             */  
long i,j,l,f,id,m,k,row,flag1,flag2;                /* counters                   */
double sum[NLBINS][FS+1];
 FILE *fp;


#define NXPIXELS 600
#define NYPIXELS 600
#define dx ((XMAX-XMIN)/(1.*NXPIXELS))
#define dy dx
 
 make_isogrid2(&isogrid);
 

 /*
   printf("%.4e %.4e\n",c[10].x,c[10].y);
   printf("rmax: %.4e NXPIXELS: %.4e XMIN %.4e\n",R_MAX,1.*NXPIXELS,1.*XMIN);
   printf("%.4e %.4e\n",isogrid.x[50][50],isogrid.y[50][50]); 
 */
 

 /**** write isodens map (ready for IDL) ****/
  
 printf("%.4e %.4e\n",1.*R_MAX_CORE,XMIN);      
 
 if (( fp = fopen("isodens.dat", "w")) == NULL)
   {fprintf(stderr, "FAILED TO OPEN isodens.dat FILE.\n");exit(0);}	       
 
 for (i=NYPIXELS-1 ; i>=0; i--)
   {
     for (j =0; j<NXPIXELS; j++) 
       {		       		  
	 r=isogrid.r[j][i];
	 /*  printf("%.4e %.4e ",1.*R_MAX_CORE/AU,XMIN/AU);printf("r=%.4e\n",r/AU);*/
	 if (isogrid.y[j][i]==0) theta=M_PI/2.;
	 else
	   {
	     if (isogrid.y[j][i]>0)
	       theta=atan(fabs(isogrid.x[j][i]/isogrid.y[j][i]));
	     else
	       theta=M_PI-atan(fabs(isogrid.x[j][i]/isogrid.y[j][i]));
	   }	 	
	   
	 if (r>R_MAX_CORE) fprintf(fp," %8.5e ",0.);

	 else
	   fprintf(fp," %8.5e  ",log10((1/(MU*M_PROTON))*dens(r,theta,0)));	 
	 
       }
     
   }
 
 fprintf(fp,"\n");
	       
fclose(fp);




return 0;

}


/******* make_isogrid *******************************************************/

void make_isogrid2(isogrid_ptr)
     
     struct isogrid_type2 *isogrid_ptr;
     
      
{
  int i,j,f,id,l,m;
  
/* printf("dx %.3e   dy %.3e\n",dx,dy); */
  

  for (i=0;i<NXPIXELS;i++)
    {
      for  (j=0;j<NYPIXELS;j++)
	
     {
       (*isogrid_ptr).x[i][j]=XMIN+(i+0.5)*dx;
       (*isogrid_ptr).y[i][j]=YMAX-(j+0.5)*dy;
       (*isogrid_ptr).r[i][j]=pow(pow((*isogrid_ptr).x[i][j],2.)+pow((*isogrid_ptr).y[i][j],2.),0.5);
       /*  printf("%.4e  %.4e %.4e  \n",(*isogrid_ptr).x[i][j]/AU,(*isogrid_ptr).y[i][j]/AU,(*isogrid_ptr).r[i][j]/AU );*/
     }
    }     
 
  return;
}
