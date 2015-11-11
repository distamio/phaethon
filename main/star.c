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
#include "inout.h"

#ifndef SPH_TREE
/********* read and set star parameters manually ****************/

void set_star_parameters(star,i)
     struct star_type star[];
     int i;
{
  star[i].npackets=STAR_NPACKETS;
  star[i].output=OUTPUT_FACTOR*star[i].npackets;
  star[i].radius= STAR_RADIUS;
  star[i].temp=STAR_TEMP;
  star[i].dilution= DILUTION;


  star[i].type=MANUAL_STAR_TYPE;
  
   if (star[i].type==0) 
   {
  star[i].x=0;
  star[i].y=0;
  star[i].z=0; 
	}
else
	{
    star[i].x=STAR_X;
    star[i].y=STAR_Y;
    star[i].z=STAR_Z;
	}
	
  star[i].mass=STAR_MASS; 
  star[i].sedfile=STAR_SED_FILE;
 /* computed in the program */

  star[i].R=sqrt(pow(star[i].x,2)+pow(star[i].y,2)+pow(star[i].z,2));

 
 if (star[i].type==3 && star[i].R<R_MAX) printf("\n:: WARNING :: Your star is inside the cloud !!!\n");	  
 	  

  return;
}

#endif


#ifdef SPH_TREE
#ifndef AUTO
/********* read and set star parameters manually ****************/

void set_star_parameters(star,i)
     struct star_type star[];
     int i;
{
  star[i].npackets=STAR_NPACKETS;
  star[i].output=OUTPUT_FACTOR*star[i].npackets;
  star[i].radius= STAR_RADIUS;
  star[i].temp=STAR_TEMP;
  star[i].dilution= DILUTION;

  star[i].rcells=STAR_RCELLS;
  star[i].lcells=STAR_LCELLS;
  star[i].mcells=STAR_MCELLS;

#ifdef STAR_DISC
  star[i].sg=STAR_GRID;
  star[i].disc=STAR_DISC;
  
  star[i].disc_mass=SMARTIE_MASS;
  star[i].disc_z=SMARTIE_z;
  star[i].disc_R= SMARTIE_R ;
  star[i].disc_kx=SMARTIE_kx ;
  star[i].disc_ky=SMARTIE_ky;
  star[i].disc_kz=SMARTIE_kz ;

#endif
    star[i].x=STAR_X;
    star[i].y=STAR_Y ;
    star[i].z=STAR_Z ; 
    star[i].mass=STAR_MASS; 
  

  /* computed in the program */
  star[i].R=sqrt(pow(star[i].x,2)+pow(star[i].y,2)+pow(star[i].z,2));
  star[i].r_dust=((star[i].radius/2.)*pow(star[i].temp/DUST_TEMP,3));
  star[i].lum=(star[i].dilution*4*M_PI*pow(star[i].radius,2)*sigma_CONST*pow(star[i].temp,4));
  star[i].disc_alpha=star[i].disc_z/star[i].disc_R;
  star[i].disc_theta_opening=F*star[i].disc_alpha;
  printf(":: Star lum   :%.5f Lsun\n",star[i].lum/SOLAR_LUM);
  printf(":: Star.Rdust :%.4e (AU)\n",star[i].r_dust/AU);

  star[i].r_star_grid=0;
  star[i].firstcellid=0;
  star[i].lastcellid=0;
  star[i].tint=5000;
  star[i].dint=100;
  return;
}
#endif
#endif

#ifdef SPH_TREE
#ifdef AUTO
/********* read and set star parameters from SinkInfo file ****************/

void set_star_parameters_auto(star,i,dragon_step,dragon_time,c,ph)
     struct star_type star[];
     int i;
     long dragon_step;
     double dragon_time;
     struct cell_type c[];
     struct photon_packet_type ph[];
{
  
  char buff[1000];
  int p;
  short int flag=0;
  double t1=0.,t2=0.,m1=0.,m2=0.;
  FILE *fp;

   struct sink_type
   {
     int    id;
     int    index;
     double x;
     double y;
     double z; 
     double vx;
     double vy;
     double vz;
     double hx;
     double hy;
     double hz;
     double mass;
     double time;
     int step;
   };

   struct sink_type sink;
   
 printf ("Reading sink data...\n");

   if (( fp= fopen("SinkInfo.data", "r")) == NULL)
     {
       fprintf(stderr, "FAILED TO OPEN SinkInfo.data FILE. \n");
       exit(0);
     }
   
#define MAX_NUMBER_OF_SPH_STEPS 1000000    
#define MAX_NUMBER_OF_SPH_STEPS 1000000    


   for (p = 0; (p <MAX_NUMBER_OF_SPH_STEPS) & getdata(fp, buff); p=p+1)
     {              
       sscanf(buff,"%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
              &(sink.time),
              &(sink.step),
              &(sink.id),
              &(sink.index),
              &(sink.mass),
              &(sink.x),
              &(sink.y),
              &(sink.z),
              &(sink.vx),
              &(sink.vy),
              &(sink.vz),
              &(sink.hx),
              &(sink.hy),
              &(sink.hz));


       if ((sink.id==i+1) && (sink.step>dragon_step)) {flag=1;break;}
     }

   fclose(fp);





  if (flag==0) error("PROBLEM IN READING SinkInfo.data!!!\n\n");
  
  star[i].sphtime=sink.time;
  star[i].sphstep=sink.step;
  star[i].sphid=sink.id;
  star[i].sphindex=sink.index;
  star[i].mass=sink.mass*SOLAR_MASS;
  /*  star[i].x=sink.x*AU;
      star[i].y=sink.y*AU;
      star[i].z=sink.z*AU; */
  star[i].vx=sink.vx;
  star[i].vy=sink.vy;
  star[i].vz=sink.vz;
  star[i].disc_kx=sink.hx/sqrt(sink.hx*sink.hx+sink.hy*sink.hy+sink.hz*sink.hz);
  star[i].disc_ky=sink.hy/sqrt(sink.hx*sink.hx+sink.hy*sink.hy+sink.hz*sink.hz);
  star[i].disc_kz=1-sqrt(star[i].disc_kx*star[i].disc_kx+star[i].disc_ky*star[i].disc_ky);

/* find accretion rate */
 
 if (( fp= fopen("SinkInfo.data", "r")) == NULL)
     {
       fprintf(stderr, "FAILED TO OPEN Sink.Info.data FILE. \n");
       exit(0);
     }
  
  for (p = 0; (p <MAX_NUMBER_OF_SPH_STEPS) & getdata(fp, buff); p=p+1)
     {              
       sscanf(buff,"%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
              &(sink.time),
              &(sink.step),
              &(sink.id),
              &(sink.index),
              &(sink.mass),
              &(sink.x),
              &(sink.y),
              &(sink.z),
              &(sink.vx),
              &(sink.vy),
              &(sink.vz),
              &(sink.hx),
              &(sink.hy),
              &(sink.hz));

       /* error ("look at dens.c!\n");*/
       if ((dragon_time-sink.time)>= 0.000004) {m1=sink.mass; t1=sink.time;}
       if ((dragon_time-sink.time)>=-0.0004) {m2=sink.mass; t2=sink.time;} /* 0.0004 */
       /* if ((dragon_time-sink.time)>= 0.0000004) {m1=sink.mass; t1=sink.time;}
          if ((dragon_time-sink.time)>=-0.0004) {m2=sink.mass; t2=sink.time;} fst */
     }

   fclose(fp);

   printf("t1:%.4e (Myr) m1:%.4e(Msun) \nt2:%.4e (Myr) m2:%.4e (Msun)\n",t1,m1,t2,m2);
   star[i].mdot=(m2-m1)/(t2-t1);

   star[i].radius= STAR_RADIUS;

   /* compute star temperature */
    
   star[i].lum=SOLAR_LUM*pow(star[i].mass/SOLAR_MASS,3)+(G_CONST*star[i].mass/star[i].radius*star[i].mdot*SOLAR_MASS/MYR);

     
   star[i].temp=pow(((SOLAR_LUM*pow(star[i].mass/SOLAR_MASS,3)+
                     (G_CONST*star[i].mass/star[i].radius)*star[i].mdot*SOLAR_MASS/MYR)/
                    (4*M_PI*star[i].radius*star[i].radius*sigma_CONST)),0.25);


  star[i].temp=STAR_TEMP;
  star[i].lum=4*M_PI*pow(star[i].radius,2.)*sigma_CONST*pow(star[i].temp,4.);
  printf("Setting star temp to %.2f\n", star[i].temp);
  printf("Setting star lum to %.2f\n", star[i].lum/SOLAR_LUM);

  star[i].npackets=STAR_NPACKETS;
star[i].output=OUTPUT_FACTOR*star[i].npackets;
  star[i].dilution= DILUTION;
  
  star[i].rcells=STAR_RCELLS;
  star[i].lcells=STAR_LCELLS;
  star[i].mcells=STAR_MCELLS;
  
  star[i].sg=STAR_GRID;
  star[i].disc=STAR_DISC;



  star[i].cell=find_star_cell(ph,0,c,star,i);
  
  printf("%d star's cell        :%4d\n",i,star[i].cell);
  printf("%d star's cell density:%.4e\n",i,c[star[i].cell].dens); 
  printf("%d star's cell size   :%.4e \n",i,c[star[i].cell].size);


  star[i].disc_alpha=DISC_ALPHA;
  star[i].disc_theta_opening=F*star[i].disc_alpha;
  star[i].disc_R=0.5*STAR_GRID_FACTOR*c[star[i].cell].size;
  star[i].disc_z=star[i].disc_alpha*star[i].disc_R;
  star[i].R=sqrt(pow(star[i].x,2)+pow(star[i].y,2)+pow(star[i].z,2));
  star[i].r_dust=((star[i].radius/2.)*pow(star[i].temp/DUST_TEMP,3));


  star[i].disc_mass=c[star[i].cell].dens*(4/3.)*M_PI*pow(STAR_GRID_FACTOR*star[i].disc_R,3);
  
  star[i].r_star_grid=0;
  star[i].firstcellid=0;
  star[i].lastcellid=0;
  star[i].tint=5000;
  star[i].dint=100;
  
  printf("\n---------- STAR %d  -----------\n",i);
  printf("Star mass: %.4e (Msun)\n",star[i].mass/SOLAR_MASS);
  printf("Star x   :% .4e (AU)\n",star[i].x/AU   );
  printf("Star y   :% .4e (AU)\n",star[i].y/AU   );
  printf("Star z   :% .4e (AU)\n",star[i].z/AU   );
  printf("Star x   :% .15e (cm)\n",star[i].x  );
  printf("Star y   :% .15e (cm)\n",star[i].y  );
  printf("Star z   :% .15e (cm)\n",star[i].z  );
  printf("Star kx  :% .15e \n",star[i].disc_kx  );
  printf("Star ky  :% .15e \n",star[i].disc_ky  );
  printf("Star kz  :% .15e \n",star[i].disc_kz  );
  printf("Star step: %d \n",star[i].sphstep);
  printf("Star Mdot:% .4e (Msun/Myr)\n",star[i].mdot);
  printf("Star Temp: %.0f (K)\n",star[i].temp);


  printf(":: Star lum   :%.5f Lsun\n",star[i].lum/SOLAR_LUM);
  printf(":: Star.Rdust :%.4e (AU)\n",star[i].r_dust/AU);
return;
}

#endif
#endif

/**** set ISRF parameters *********************************************************************/
void set_isrf_parameters(star,i)
     struct star_type star[];
     int i;
 
{
  
  star[i].npackets=ISRF_NPACKETS;
  star[i].radius= R_ISRF;
  star[i].output=OUTPUT_FACTOR*star[i].npackets;
  star[i].tint=50;
  star[i].dint=1;
  
  return;
}

/**** set STAR ISRF parameters ******************************************************************/
void set_star_isrf_parameters(star,i)
     struct star_type star[];
     int i;
 
{
  
  star[i].npackets=ISRF_NPACKETS;
  star[i].radius= R_ISRF;
  star[i].output=OUTPUT_FACTOR*star[i].npackets;
  star[i].tint=50;
  star[i].dint=1;
  star[i].temp=STAR_ISRF_TEMP;
  star[i].dilution= DILUTION;
  star[i].lum=(star[i].dilution*4*M_PI*pow(star[i].radius,2)*sigma_CONST*pow(star[i].temp,4));
  return;
}




