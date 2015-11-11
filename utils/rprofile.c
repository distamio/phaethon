/* 

calculate radial profiles at a specific x pixel (i.e. vertical cuts in the isophotal maps)

       
output files: radial.$lambda.$id.$l.$m.$pixel

needs to be compil and run locally at the directory under study

sm file to use: rprof.sm

compile: gcc -g -o rprof  rprofile.c /home/athena/spxstd/npha/nphaethon/main/inout.c -lm

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>
#include "constants.h"
#include "params.h"
#include "params.grid.h"
#include "params.opa.h"

#include "structs.h"
#include "inout.h"
#include "phaethon.functions.h"
#include "photons.h"
#include "spec.functions.h"
#include "functions.h"
#include "phaethon.h"


#define FS1 7
#define FS2 7

char      filename[20],filename1[20],filename2[20],filename3[30];

int main()

{
double photon_i=0,dlambda;         /* photon counter             */  
int i,j,l,f,id,m,k,row,col;                /* counters                   */
double sum[NLBINS][NMBINS][FS+1];
double lambda_obs[FS]; /* wavelengths to construct isophotal maps */
int iNLBINS,iNMBINS;

FILE *fp;

iNLBINS=NLBINS;
iNMBINS=NMBINS;

#ifdef SPHERICAL
iNLBINS=1;
iNMBINS=1;
#endif

#ifdef DUST_MIX
  read_grid_params (opacity,opacity2);
#else
  read_grid_params (opacity);
#endif


 // read run parameters (uses info from the grid) 

    read_params ();


 // make grid



 for (i=0;i<FS;i++)
     {
         switch (i)
           {
           case 0:lambda_obs[0]=F0;break;
           case 1:lambda_obs[1]=F1;break;
           case 2:lambda_obs[2]=F2;break;
           case 3:lambda_obs[3]=F3;break;
           case 4:lambda_obs[4]=F4;break;
           case 5:lambda_obs[5]=F5;break;
           case 6:lambda_obs[6]=F6;break;
           case 7:lambda_obs[7]=F7;break;
           case 8:lambda_obs[8]=F8;break;
           case 9:lambda_obs[9]=F9;break;
           }

       }

printf("NLBINS: %d NMBINS: %d\n",iNLBINS,iNMBINS);

/* reads isogrids files */

 for (l=0;l<iNLBINS;l++)
   for (m=0;m<iNMBINS;m++)
     {
       for (k=FS1;k<=FS2;k++)
	 {
	   sum[l][m][k]=0.;
	   read_isogrid(&isogrid,l,m);

           for (i=0; i<NXPIXELS; i++)
	     {
	       for (j=NYPIXELS-1; j>=0; j--) 
		 {
		   sum[l][m][k]=sum[l][m][k]+isogrid.intensity[i][j][k][2][l][m];
		  
		 }
	     }	   
	 }
     }


 printf("\n Integrated intensity (MJy/sr) at different viewing angles (just from the source)\n(lambda, l, m, intensity) \n\n");
 for (k=FS1;k<=FS2;k++)
   {   
     for (l=0;l<iNLBINS;l++)
       for (m=0;m<iNMBINS;m++)
       {
	 printf(" %d %d %d sum=%.3e\n",k,l,m,sum[l][m][k]);
       }
       printf("----------------------\n");
   }



 for (l=0;l<iNLBINS;l++)
     for (m=0;m<iNMBINS;m++)
       { 
 	for (row=63;row<64+0*NYPIXELS;row=row+1)
   	{     
     	 for (id=0;id<=4;id++)
       		{


	 	for (k=FS1;k<=FS2;k++)
	   		{

  	     		sprintf(filename3,"%s.%d.%d.%d.%d.%d.dat","radial",k,id,l,m,row);


	     		if (( fp = fopen(filename3, "w")) == NULL)
	       		{printf("FAILED TO OPEN %s FILE. \n",filename3);}

	     		for (i=NYPIXELS-1 ; i>=0; i--)
	      			 {

		 		fprintf(fp," %d % .3e",NYPIXELS-1-i,(isogrid.intensity[row-3][i][k][id][l][m]+isogrid.intensity[row-2][i][k][id][l][m]+isogrid.intensity[row-1][i][k][id][l][m]+isogrid.intensity[row][i][k][id][l][m]+isogrid.intensity[row+1][i][k][id][l][m]+isogrid.intensity[row+2][i][k][id][l][m]+isogrid.intensity[row+3][i][k][id][l][m])/7.);
		 
		 		fprintf(fp,"\n");
	       			}	
	     
	     		fclose(fp);
	   		}	
       		}       
   	}
       }



 return 0;
 
}
