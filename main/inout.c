#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "inout.h" 
#include "constants.h"
#include "params.h"
#include "phaethon.params.h"

#define MAX_LINE     500

/***** function read_photons  *****************************************/

/* warning: this does not change star params. i need ptr to do this   */

/* void read_photons(p)
    struct photon_packet_type p[];
{
   int i;
   char buff[MAX_LINE];
   FILE *fp;

   if (( fp = fopen("inphotons.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN photons FILE. \n");
         exit(0);
         }

for (i = 0; (i <NPACKETS) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d",
              &(p[i].x),
              &(p[i].y),
              &(p[i].z),
              &(p[i].kx),
              &(p[i].ky),
              &(p[i].kz),
              &(p[i].f),
              &(p[i].lum),
              &(p[i].tau),
	      &(p[i].starlum), 
	      &(p[i].abs),
	      &(p[i].scat),
	      &(p[i].lcell),
	      &(p[i].id)
              );       
  }

   fclose(fp); 
  
   return;
}
*/
/***** function write_photons *****************************************/
/*
void write_photons(p)

   struct photon_packet_type p[];

      
{
   int i;
   FILE *fp;
 

   
   if (( fp = fopen("inphotons.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN inphotons.dat FILE. \n");
         exit(0);
         }


fprintf(fp, "@ x,y,z,kx,ky,kz,f,lum,tau,starlum,abs,scat,lcell,id \n");


for (i=0;i<NPACKETS;i++)

fprintf(fp, "% .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e %.3f %5d %5d %5d %2d\n", 
              p[i].x,
              p[i].y,
              p[i].z,
              p[i].kx,
              p[i].ky,
              p[i].kz,
              p[i].f,
              p[i].lum,
              p[i].tau,
	      p[i].starlum,
	      p[i].abs,
              p[i].scat,
              p[i].lcell,
              p[i].id
               );               
     

   fclose(fp);
  
   return;
}
*/
/**************** define getline **************************************/

int mygetline(fp, buff)
  FILE *fp;
  char buff[];
{
   int i = 0;
   char c;

   while (TRUE)
    {
      c = getc(fp);
      if (feof(fp)) 
         {
           return(FALSE);
         }	  
         else
         {
             buff[i] = c;	   
	     if (c == '\n')
               {
                buff[i] = '\0';
                return(TRUE);
     	       }	
              i++;
	  }
    }

}

/************** define getdata *****************************************/

getdata(fp, buff)
 FILE *fp;
 char buff[];
{
  
  while (mygetline(fp, buff))
  {
     if (buff[0] != '@' && buff[0] != 'l') 
     {
	 return(TRUE);
     }
  }
  return(FALSE);
}



/***** function read_treebo  *****************************************/
/*
void read_treebo(p,nparticles_ptr)

    struct particle_type p[];
    int  *nparticles_ptr;

{
   int c;
   char buff[MAX_LINE];
   FILE *fp;

  
   if (( fp = fopen("treebo", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN treebo FILE. \n");
         exit(0);
         }

   for (c = 0; (c <MAX_PARTICLE) & getdata(fp, buff); c++)
     {
       sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     &(p[c].x),
                     &(p[c].y),
                     &(p[c].z),
                     &(p[c].vx),
                     &(p[c].vy),
                     &(p[c].vz),
                     &(p[c].rho),
                     &(p[c].mass),
                     &(p[c].temp),
                     &(p[c].iwas),
                     &(p[c].h),
                     &(p[c].lum)
                     );
       *nparticles_ptr=(*nparticles_ptr)+1;
      }

   fclose(fp);

   return;
}

*/
/****************** Define write_treebo *********************************/


void write_treebo(p,nparticles)

    struct particle_type p[];
    long int nparticles;
{
   int c;
   FILE *fp;

   fp = fopen("treebo.out", "w");

   for (c = 0; c <nparticles; c++)
     {
     fprintf(fp, 
 "% .6e % .6e % .6e % .6e % .6e % .6e % .6e % .6e % .6e %d % .6e % .1e\n", 
                     p[c].x,                
                     p[c].y,
                     p[c].z,
                     p[c].vx,
                     p[c].vy,
                     p[c].vz,
                     p[c].rho,
                     p[c].mass,
                     p[c].temp,
                     p[c].iwas,
                     p[c].h,
                     p[c].lum
                 );
      }
   fclose(fp);
   return;
}


/****************** write cells*********************************/

void write_cells(c, ncells)

 struct cell_type c[];
 int ncells;

  {

   int i;
   FILE *fp;
 

   /* Open tfile */
   if (( fp = fopen("incells.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN cells.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "@                   x                         y                        z            nabs nscat    abslum      temp    dens        mass       cell_tau  next_n\n");
 
for (i=0;i<ncells;i++)

fprintf(fp, " % .20e % .20e % .20e %4ld %4ld % .4e %6.1f % .4e % .4e % .4e %4ld\n", 

      c[i].x,
      c[i].y, 
      c[i].z,     
      c[i].nabs,
      c[i].nscat,
      c[i].abslum,
      c[i].temp,
      c[i].dens,
      c[i].mass,c[i].tau, 
      c[i].next_n );               
     
   fclose(fp);
  
   return;
}


/****************** read cells *********************************/

void read_cells(c, ncells_ptr)

 struct cell_type c[];
 long *ncells_ptr;

  {
   char buff[MAX_LINE];
   int i;
 
   FILE *fp;
 

   /* Open tfile */
   if (( fp = fopen("incells.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN incells.dat FILE. \n");
         exit(0);
         }

for (i = 0; (i <MAX_CELLS) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%lf %lf %lf %ld %ld %lf %lf %lf %lf %lf %ld",
          &(c[i].x),
          &(c[i].y),
          &(c[i].z), 
          &(c[i].nabs),
          &(c[i].nscat),
          &(c[i].abslum),
          &(c[i].temp),
          &(c[i].dens),
          &(c[i].mass),
          &(c[i].tau),
          &(c[i].next_n) );    
               
     *ncells_ptr=(*ncells_ptr)+1;
      }

   fclose(fp);
  

   /* routine to find A_V

   double tau[NRCELLS+1];
   FILE *fp;
   if (( fp = fopen("Av.r.dat", "r")) == NULL)
   {
   fprintf(stderr, "FAILED TO OPEN incells.dat FILE. \n");
   exit(0);
   }
   
   tau[0]=0;
   
   for (i = 1; i<=NRCELLS; i++) 
   
   tau[i]=c[NRCELLS-i-1].tau+tau[i-1];
   
   
   for (i =0; i <=NRCELLS; i++) 
   
   printf("A_v=%.3f  r=%.4e \n", tau[i]*1.086, c[NRCELLS-i-1].x/AU);
   
   exit(0);
   */

   return;

}


/****************** read cells *********************************/

void read_cells_xyz(c, ncells_ptr)

 struct cell_type c[];
 long *ncells_ptr;

  {
   char buff[MAX_LINE];
   int i;
   double junk;

   FILE *fp;
 

   /* Open tfile */
   if (( fp = fopen("incells.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN incells.dat FILE. \n");
         exit(0);
         }

for (i = 0; (i <MAX_CELLS) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
          &(c[i].x),
          &(c[i].y),
          &(c[i].z), 
          &junk,
          &junk,
          &junk,
          &junk,
          &junk,
          &junk,
          &junk, &junk);         
     *ncells_ptr=(*ncells_ptr)+1;
      }

   fclose(fp);
  

   return;

}


/****************** write grid*********************************/

void write_grid(grid)

 struct grid_type grid;

  {

   FILE *fp;
 

   /* Open tfile */
   if (( fp = fopen("grid.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN grid.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "@dx,dy,dz\n");
fprintf(fp, "% .3e % .3e % .3e\n", 
grid.dx,grid.dy,grid.dz);               
     
   fclose(fp);
  
   return;
}


/***** function read_grid  *****************************************/

void read_grid(grid_ptr)

     struct grid_type *grid_ptr;

{
   int c;
   char buff[MAX_LINE];
   FILE *fp;

   /* Open the treebo file */
   if (( fp = fopen("grid.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN grid.dat FILE. \n");
         exit(0);
         }

   for (c = 0; (c <2) & getdata(fp, buff); c++)
     {
       sscanf(buff, "%lf %lf %lf",
              &(*grid_ptr).dx,&(*grid_ptr).dy,&(*grid_ptr).dz);      
      }

   fclose(fp);

   return;
}


/***** function read_opacity  *****************************************/

/*  opacity.dat file should contain lambda, absorption, scattering    */

/*  lambda should be read in microns, rest in cgs (cm^2/gr)           */

void read_opacity(opacity_ptr)

    struct opacity_type *opacity_ptr;

{
   int i;
   FILE *fp;

   /* Open  file */
   if (( fp = fopen("opacity.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
         exit(0);
         }

for (i = 0;i <OPADIM; i++)
  {    
   fscanf(fp,"%lf %lf %lf\n",&((*opacity_ptr).lambda[i]),
           &((*opacity_ptr).abs[i]),& (*opacity_ptr).scat[i]);

   (*opacity_ptr).scat[i]=SCATTER_MULTIPLY*((*opacity_ptr).scat[i]);
   (*opacity_ptr).abs[i]=((*opacity_ptr).abs[i]);
   
   (*opacity_ptr).albedo[i]=
     ((*opacity_ptr).scat[i])/((*opacity_ptr).scat[i]+ (*opacity_ptr).abs[i]);
  }  
     
   fclose(fp);

   return;
}

#ifdef DUST_MIX

/***** function read_opacity2  *****************************************/

/*  opacity.dat file should contain lambda, absorption, scattering    */

/*  lambda should be read in microns, rest in cgs (cm^2/gr)           */

void read_opacity2(opacity_ptr)

    struct opacity_type *opacity_ptr;

{
   int i;
    char buff[MAX_LINE]; 
  
   FILE *fp;

   /* Open  file */
   if (( fp = fopen("opacity2.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity2.dat FILE. \n");
         exit(0);
         }

for (i = 0;i <OPADIM; i++)
  {    
   fscanf(fp,"%lf %lf %lf\n",&((*opacity_ptr).lambda[i]),
           &((*opacity_ptr).abs[i]),& (*opacity_ptr).scat[i]);
  
   (*opacity_ptr).scat[i]=SCATTER_MULTIPLY*((*opacity_ptr).scat[i]);
   (*opacity_ptr).abs[i]=((*opacity_ptr).abs[i]);

   (*opacity_ptr).albedo[i]=
     ((*opacity_ptr).scat[i])/((*opacity_ptr).scat[i]+ (*opacity_ptr).abs[i]);
  }  
     
   fclose(fp);

   return;
}

#endif

/***** function write_opacity  *****************************************/


void write_opacity(opacity)

    struct opacity_type opacity;

{
   int i;

   FILE *fp;

   /* Open  file */
   if (( fp = fopen("opacity.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN opacity.dat FILE. \n");
         exit(0);
         }

for (i = 0;i <OPADIM; i++)
     
      fprintf(fp,"%.5e %4e %4e\n",opacity.lambda[i],opacity.abs[i],opacity.scat[i]);
      
   fclose(fp);

   return;
}

/***** function record_photons *****************************************/

/* 
void record_photons(p)

   struct photon_packet_type p[];
      
{
   int i;
   FILE *fp;
 
   
   if (( fp = fopen("outphotons.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN out_photons.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "@ x,y,z,kx,ky,kz,f,lum,tau,starlum,abs,scat,lcell,id \n");

for (i=0;i<NPACKETS;i++)

fprintf(fp, "% .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e %.4e %5d %5d %5d %2d\n", 
              p[i].x,
              p[i].y,
              p[i].z,
              p[i].kx,
              p[i].ky,
              p[i].kz,
              p[i].f,
              p[i].lum,
              p[i].tau,
	      p[i].starlum,
	p[i].abs,p[i].scat,p[i].lcell,p[i].id
               );                  

   fclose(fp);
  
   return;
}
*/
/****************** record cells*********************************/

void record_cells(c, ncells)

 struct cell_type c[];
 int ncells;

  {
   int i;
   FILE *fp;
 
   /* Open tfile */
   if (( fp = fopen("outcells.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN outcells.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "@    x                 y               z           nabs nscat  abslum      temp    dens       mass       cell_tau  next_n\n"); 
for (i=0;i<ncells;i++)

fprintf(fp, " % .7e % .7e % .7e %4ld %4ld % .4e %6.2f % .4e % .4e % .4e  %4ld\n", 

	c[i].x,
	c[i].y,
	c[i].z, 
	c[i].nabs,
	c[i].nscat,
	c[i].abslum,
	c[i].temp,
	c[i].dens,
	c[i].mass,
	c[i].tau, c[i].next_n);               
     
 fclose(fp);
 
 return;
  }


#ifdef SPH_TREE

/****************** record active cells*********************************/

void record_active_cells(c, ncells)

 struct cell_type c[];
 int ncells;

  {
   int i;
   FILE *fp;
 
   /* Open tfile */
   if (( fp = fopen("outcells.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN outcells.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "@    x             y             z       nabs nscat  abslum      temp    dens       mass       cell_id \n");

for (i=0;i<ncells;i++)

if (c[i].temp>0) 
  {
    fprintf(fp, " % .7e % .7e % .7e %4d %4d % .4e %6.1f % .4e % .4e %4d \n", 
	    
	    c[i].x,
	    c[i].y,
	    c[i].z, 
	    c[i].nabs,
	    c[i].nscat,
	    c[i].abslum,
	    c[i].temp,
	    c[i].dens,
	    c[i].mass,
	    i);        
       
  }
 fclose(fp); 

 return;
  }



/****************** read active cells *********************************/

void read_active_cells(c,gindex)

 struct cell_type c[];
 int gindex;

{
  char buff[MAX_LINE];
  int i;
  int index;
  double x,y,z,abslum,temp,dens,mass;
  int nabs,nscat;

  FILE *fp;
  
  
  
  if (( fp = fopen("outcells.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN outcells.dat FILE. \n");
      exit(0);
    }
  
  for (i = 0; (i <MAX_CELLS) & getdata(fp, buff); i++)
    { 
      sscanf(buff,"%lf %lf %lf %d %d %lf %lf %lf %lf %d",
	     &x,
	     &y,
	     &z,                                                  
	     &nabs,
	     &nscat,
	     &abslum,
	     &temp,
	     &dens,
	     &mass,
	     &index);       
     
      
      /* 
	 if (i==0)
	 {
	 printf("\n"); 
	 printf("  x             y            z     nabs  nscat abslum   temp    dens    mass   cell_id  nsph size\n");
	 }
	 if ((temp>8) && fabs(x)<200*AU)
	 { 
	 printf("% .4e % .4e % .4e %5d %d  %.4e %4.1f %.4e %.4e %d\n",
	 x,
	 y,
	 z,                                                  
	 nabs,
	 nscat,
	 abslum,
	 temp,
	 dens,
	 mass,
	 index);      
	 }
	 
      */

      if (index>gindex) break;
      else
	{
	  c[index].x= x;
	  c[index].y= y;
	  c[index].z= z;                                                
	  c[index].nabs= nabs;
	  c[index].nscat= nscat;
	  c[index].abslum= abslum;
	  c[index].temp= temp;
	  c[index].dens= dens;
	  c[index].mass= mass;

       }
     
     

    }
  
     
  fclose(fp);
 
  return;

}



#ifdef SPH_TREE
/****************** record active cells*********************************/

void record_potential_cell_info(c, ncells,mass)

 struct cell_type c[];
 int ncells;
 double mass;
 
{
  long i;
  double tmass=0.;
  double tnsph=0.;
  long cellcounter=0.;

  FILE *fp;
  

   for (i=0;i<ncells;i++)
     {  
       if (c[i].nsph>SPH_PARTICLES_PER_CELL)
	 {
	   tmass=tmass+c[i].mass;
	   tnsph=tnsph+c[i].nsph;
	   cellcounter=cellcounter+1;
	 }
     }
  
  /* Open tfile */
  if (( fp = fopen("potential.cells.sph.dat", "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN potential.cells.sph.dat FILE. \n");
      exit(0);
    }
  
  fprintf(fp,"@ TOTAL MASS=%.4e TOTAL MASS (real)=%.4e\n", tmass,1.*mass);
  fprintf(fp,"@ TOTAL NSPH=%.0f \n", tnsph);
  fprintf(fp,"@ Active cells: %4d\n",cellcounter);
  
  fprintf(fp, "@    x               y               z           nabs  nscat abslum   temp    dens       mass   cell_id  nsph size\n");
  
  for (i=0;i<ncells;i++)
    
    if (c[i].nsph>SPH_PARTICLES_PER_CELL)
      {
	fprintf(fp, " % .7e % .7e % .7e %4d %4d % .2e %6.1f % .4e %.4e %4d %3.0f %.2e \n", 
		
		c[i].x,	            c[i].y,	    c[i].z, 
		c[i].nabs,	    c[i].nscat,	    c[i].abslum,
		c[i].temp,	    c[i].dens,	    c[i].mass,
		i,       	    c[i].nsph,	    c[i].size);        
      }    
  

  fclose(fp);

  return;
}




/******************* record active cell info  *******************/
void record_active_cell_info(c, ncells,mass)

 struct cell_type c[];
 int ncells;
 double mass;
{
  long i;
  double tmass=0.;
  double tnsph=0.;
  long cellcounter=0.;

  FILE *fp;
  

   for (i=0;i<ncells;i++)
     {  
       if (c[i].temp>0) 
	 {
	   tmass=tmass+c[i].mass;
	   tnsph=tnsph+c[i].nsph;
	   cellcounter=cellcounter+1;
	 }
     }
  
  /* Open tfile */
  if (( fp = fopen("outcells.sph.dat", "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN outcells.sph.dat FILE. \n");
      exit(0);
    }
  
  fprintf(fp,"@ TOTAL MASS=%.4e TOTAL MASS (real)=%.4e\n", tmass,1.*mass);
  fprintf(fp,"@ TOTAL NSPH=%.0f \n", tnsph);
  fprintf(fp,"@ Active cells: %4d\n",cellcounter);
  
  fprintf(fp, "@    x               y               z           nabs  nscat abslum   temp    dens       mass   cell_id  nsph size\n");
  
  for (i=0;i<ncells;i++)
    
    if (c[i].temp>0) 
      {
	fprintf(fp, " % .7e % .7e % .7e %4d %4d % .2e %6.1f % .4e %.4e %4d %3.0f %.2e \n", 
		
		c[i].x,	            c[i].y,	    c[i].z, 
		c[i].nabs,	    c[i].nscat,	    c[i].abslum,
		c[i].temp,	    c[i].dens,	    c[i].mass,
		i,       	    c[i].nsph,	    c[i].size);        
      }    
  

  fclose(fp);

  return;
}

#endif


/******************* POSTRT OUTCELLS *******************/

#ifdef SPH_TREE

void record_active_cell_info_POSTRT(c, ncells,mass)

 struct cell_type c[];
 int ncells;
 double mass;
{
  long i;
  double tmass=0.;
  double tnsph=0.;
  long cellcounter=0.;

  FILE *fp;
  

   for (i=0;i<ncells;i++)
     {  
       if (c[i].temp>0) 
	 {
	   tmass=tmass+c[i].mass;
	   tnsph=tnsph+c[i].nsph;
	   cellcounter=cellcounter+1;
	 }
     }
  
  /* Open tfile */
  if (( fp = fopen("outcells.sph.PRT.dat", "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN outcells.sph.PRT.dat FILE. \n");
      exit(0);
    }
  
  fprintf(fp,"@TOTAL MASS=%.4e TOTAL MASS (real)=%.4e\n", tmass,1.*mass);
  fprintf(fp,"@TOTAL NSPH=%.0f \n", tnsph);
  fprintf(fp,"@Active cells: %4d\n",cellcounter);
  
  fprintf(fp, "@    x          y          z    nabs nscat abslum   temp    dens       mass   cell_id  nsph size\n");
  
  for (i=0;i<ncells;i++)
    
    if (c[i].temp>0) 
      {
	fprintf(fp, " % .7e % .7e % .7e %4d %4d % .2e %6.1f % .4e %.4e %4d %3.0f %.2e \n", 
		
		c[i].x,	            c[i].y,	    c[i].z, 
		c[i].nabs,	    c[i].nscat,	    c[i].abslum,
		c[i].temp,	    c[i].dens,	    c[i].mass,
		i,       	    c[i].nsph,	    c[i].size);        
      }    
  

  fclose(fp);

  return;
}

#endif




#endif 
/***** function read_mean_opacity  *****************************************/

void read_mean_opacity(mean_op_ptr)

    struct mean_op_type *mean_op_ptr;

{
   int i;

  
   FILE *fp;

   if (NFREQ!=FPOINTS){printf("Error: NFREQ must be equal to FPOINTS\n");exit(0);}
   /* Open  file */
   if (( fp = fopen("mopacity.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN mopacity.dat FILE. \n");
         exit(0);
         }

   for (i = 0;i <TEMPDIM; i++)
   
     fscanf(fp,"%lf %lf\n", & ((*mean_op_ptr).temp[i]),& ((*mean_op_ptr).opacity[i])); 
     
   fclose(fp);

   return;
}

#ifdef DUST_MIX
/***** function read_mean_opacity2  *****************************************/

void read_mean_opacity2(mean_op_ptr)

    struct mean_op_type *mean_op_ptr;

{
   int i;
    char buff[MAX_LINE]; 
  
   FILE *fp;

   /* Open  file */
   if (( fp = fopen("mopacity2.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN mopacity2.dat FILE. \n");
         exit(0);
         }

   for (i = 0;i <TEMPDIM; i++)
   
     fscanf(fp,"%lf %lf\n", & ((*mean_op_ptr).temp[i]),& ((*mean_op_ptr).opacity[i])); 
     
   fclose(fp);

   return;
}

#endif

/***** function write_mean_opacity  *****************************************/

/* cgs units */
void write_mean_opacity(mean_op)

    struct mean_op_type mean_op;

{
   int i;
  
   FILE *fp;

   /* Open  file */
   if (( fp = fopen("mopacity.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN mopacity.dat FILE. \n");
         exit(0);
         }
  

   for (i =0;i<TEMPDIM; i++)
   
   fprintf(fp,"%6.2f %.4e\n", mean_op.temp[i],mean_op.opacity[i]); 
     
   fclose(fp);

   return;
}
/***** function read_outphotons  *****************************************/
/*
void read_outphotons(p)
    struct photon_packet_type p[];
{
   int i;
   char buff[MAX_LINE];
   FILE *fp;

   if (( fp = fopen("outphotons.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN photons FILE. \n");
         exit(0);
         }

for (i = 0; (i <NPACKETS  ) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d",
              &(p[i].x),
              &(p[i].y),
              &(p[i].z),
              &(p[i].kx),
              &(p[i].ky),
              &(p[i].kz),
              &(p[i].f),
              &(p[i].lum),
              &(p[i].tau),
	      &(p[i].starlum),
	      &(p[i].abs),
	      &(p[i].scat),
	      &(p[i].lcell),
	      &(p[i].id)

              );       
  }

   fclose(fp);
  
   return;
}

*/
/****************** read outcells *********************************/

void read_outcells(c, ncells_ptr)

 struct cell_type c[];
 long *ncells_ptr;

{
  char buff[MAX_LINE];
  int i;
  double junk;
  FILE *fp;
  
  
  
  if (( fp = fopen("outcells.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN outcells.dat FILE. \n");
      exit(0);
    }
  
  for (i = 0; (i <MAX_CELLS) & getdata(fp, buff); i++)
    { 
      sscanf(buff,"%lf %lf %lf %ld %ld %lf %lf %lf %lf %lf %ld",
	     &junk,
	     &junk,
	     &junk,                                                  
	     &(c[i].nabs),
	     &(c[i].nscat),
	     &(c[i].abslum),
	     &(c[i].temp),
	     &(c[i].dens),
	     &(c[i].mass),
	     &(c[i].tau), &c[i].next_n);         
      *ncells_ptr=(*ncells_ptr)+1;
    }
  
    
  fclose(fp);


  return;

}



/***** function record_sed_info*****************************************/

/*
void record_sed_info(p)

   struct photon_packet_type p[];
      
{
   int i;
   FILE *fp;
 
  
   if (( fp = fopen("outsed.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN outsed.dat FILE. \n");
         exit(0);
         }

fprintf(fp, "fbin l  m id lcell\n");

for (i=0;i<NPACKETS;i++)

fprintf(fp,"%3d %2d %2d %2d %3d\n", 
              p[i].n,
              p[i].l,
              p[i].m,
	      p[i].id,
              p[i].lcell          
               );                  

   fclose(fp);
  
   return;
}
*/
/***** function read_sed_info *****************************************/
/*
void read_sed_info(p)
    struct photon_packet_type p[];
{
   int i;
   char buff[MAX_LINE];
  
   FILE *fp;

   if (( fp = fopen("outsed.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN outsed FILE. \n");
         exit(0);
         }

for (i = 0; (i <NPACKETS  ) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%d %d %d %d %d",
              &(p[i].n),
              &(p[i].l),
              &(p[i].m),
	      &(p[i].id),
              &(p[i].lcell)
             
              );       
  }

   fclose(fp);
  
   return;
}
*/
/***** function write_ftable *****************************************/

void write_ftable(reemit)

     struct reemit_type reemit;

{
   int k;
   FILE *fp;
 
   /* Open tfile */
   if (( fp = fopen("ftable.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ftable.dat FILE. \n");
         exit(0);
         }

  /* write ftable.dat (the reemission frequencies) */

for (k=0; k<NFREQ ; k++)  fprintf(fp,"%.4e\n", reemit.freq[k]);         
	   
fclose(fp);

   return;
}



/***** function read_ftable *****************************************/

void read_ftable(reemit_ptr)

     struct reemit_type *reemit_ptr;
{
   int i;
   FILE *fp;

   if (( fp = fopen("ftable.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ftable.dat FILE. \n");
         exit(0);
         }

for (i = 0; i <NFREQ; i++)
  { 
  fscanf(fp,"%lf", &(*reemit_ptr).freq[i]);       
  }

   fclose(fp);
  
   return;
}

#ifdef DUST_MIX
/***** function read_ftable2 *****************************************/

void read_ftable2(reemit2_ptr)

     struct reemit_type *reemit2_ptr;
{
   int i;
   char buff[MAX_LINE];
  
   FILE *fp;

   if (( fp = fopen("ftable2.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ftable.dat FILE. \n");
         exit(0);
         }

for (i = 0; i <NFREQ; i++)
  { 
  fscanf(fp,"%lf", &(*reemit2_ptr).freq[i]);       
  }

   fclose(fp);
  
   return;
}

#endif

/****************** read ISRF*********************************/

/* read a table where lambda decreases and put it in a table
   where lambda increases */

void read_ISRF(bisrf_ptr)

     struct btype *bisrf_ptr;

  {
   char buff[MAX_LINE];
   int i;
   FILE *fp;


   /* Open tfile */
   if (( fp = fopen("isrf.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN isrf.dat FILE. \n");
         exit(0);
         }
   (*bisrf_ptr).dim=0;

   for (i = 0; (i <ISRF_MAX_DIM) & (getdata(fp, buff)); i++)
     {
       sscanf(buff,"%lf %lf",&(*bisrf_ptr).lambda[i],&(*bisrf_ptr).intensity[i]);        
       (*bisrf_ptr).dim=(*bisrf_ptr).dim+1;

  if ((*bisrf_ptr).lambda[i]<500) (*bisrf_ptr).intensity[i]=(*bisrf_ptr).intensity[i]*ISRF_MULTIPLY;
     }

   fclose(fp);
 
   return;
}

/****************** read STAR SED*********************************/

void read_star_sed(star,s)
     struct star_type star[];
     int s;
  {
    char buff[MAX_LINE];
    int i;
    char filename[40];
    double junk,junk2,help1[ISRF_MAX_DIM],help2[ISRF_MAX_DIM];
    FILE *fp;
    
    sprintf(filename,"star%d.sed.dat",s);
    
    /* Open tfile */
   if (( fp = fopen(filename, "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN star%d.sed.dat FILE.\n",s);
         exit(0);
         }

   star[s].sed_dim=0;


// Allard COND models [lambda (in Angstrom), Flux (in ERG/CM2/S/A)]

   for (i = 0; (i <ISRF_MAX_DIM) & (getdata(fp, buff)); i++)
	
     {
       sscanf(buff,"%lf %lf",&help1[i], &help2[i]); 
       star[s].sed_dim=star[s].sed_dim+1;
     }


   fclose(fp);
  
 for (i = 0; i<star[s].sed_dim;i=i+1)
     {
      star[s].lambda[i]=help1[i]*1e-4; // transform to micron
      star[s].intensity[i]=(pow(star[s].lambda[i]*1e-4,2)/3e10)*help2[i]*1e8/M_PI; // transform flux to intensity_nu
    }

   
   return;

}

/***** function record_photons *****************************************/

void record_photons_info(p,photon_i)

   struct photon_packet_type p[];
   double  photon_i;   
{
  int i;
   FILE *fp;

   if (photon_i==0) 
     {
      fp = fopen("photons.info.dat", "w");
      fprintf(fp, "@ n/ l/ m / id \n");
      fclose(fp);
     }

   i=0;

   /* Open tfile */
   if (( fp = fopen("photons.info.dat", "a+")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN photons.info.dat FILE. \n");
         exit(0);
         }
  
fprintf(fp, "%d %d %d %d\n",               
              p[i].n,
              p[i].l,
              p[i].m,
	      p[i].id
               );                  
   
   fclose(fp);
  
   return;
}




/***** function read_outphotons_basic  *****************************************/
/*
void read_photons_info(p)
    struct photon_packet_type p[];
{
   int i;
   char buff[MAX_LINE];
   FILE *fp;

   if (( fp = fopen("bphotons.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN bphotons FILE. \n");
         exit(0);
         }

for (i = 0; (i <NPACKETS  ) & getdata(fp, buff); i++)
  { 
  sscanf(buff,"%d %d %d %d ",
              &(p[i].n),
              &(p[i].l),
              &(p[i].m),
              &(p[i].id)
              );       
  }

   fclose(fp);
  
   return;
}

*/
/***** function write_info_file ********************************/

void write_info_file(p)

   struct photon_packet_type p[];

      
{
   FILE *fp;
 

   /* Open tfile */
   if (( fp = fopen("info.dat", "w")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN info.dat FILE. \n");
         exit(0);
         }


fprintf(fp, "___________RUN INFO FILE______________\n\n");


     

   fclose(fp);
  
   return;
}




/***** function read_photons  *****************************************/

/* warning: this does not change star params. i need ptr to do this   */

/*
void read_outphotons_basic(p)
    
     struct photon_packet_type p[];
{
   long i;  
   char buff[100];
   FILE *fp;

   if (( fp = fopen("bphotons.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN bphotons FILE. \n");

         exit(0);
         }

fgets(buff,sizeof(buff),fp);

for (i = 0;i<NPACKETS; i++)
  
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %d",
 &(p[i].lx),&(p[i].ly),&(p[i].lz),&(p[i].kx),&(p[i].ky),&(p[i].kz),&(p[i].f),&(p[i].id));       
  
   fclose(fp);
  
   return;
}

*/



/****************** read isogrid *********************************/
#ifdef DO_ISOMAP
void read_isogrid(isogrid_ptr,l,m)

struct isogrid_type *isogrid_ptr;
int l,m;

{
  char junk[30];
  char buff[MAX_LINE];
  int i,id,k,j;  
  char filename[20],filename1[20],filename2[20],filename3[20];
  FILE *fp;
  



  /* read IDL files */
 
sprintf(filename,"iso");


  for (id=0;id<=4;id++) 
    { 
      for (k=0;k<FS;k++) 
	{ 
	  
	  sprintf(filename2,"%s.%d.%d.%d.%d.dat",filename,k,id,l,m);  
	  
	  if (( fp = fopen(filename2, "r")) == NULL) 
	    {
	      fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename2); 
	      exit(0); 
	    } 
	  
	  
	  for (i=NYPIXELS-1 ; i>=0; i--) 
	    { 
	      for (j =0; j<NXPIXELS; j++)  
		{ 	
		  fscanf(fp,"%lf \n",&((*isogrid_ptr).intensity[j][i][k][id][l][m]));  
		  
		} 
	    }
	   
	  fclose(fp);  
	}   
    } 
 
  return;
  
}
   
#endif

/*******************************************************************************************************************************/
void read_bbins(bbin_ptr)

struct bbin_type *bbin_ptr; 

{

  char junk[100];
  int i,id,k;  
  char filename[20],filename2[20];
  double junk2;

  FILE *fp;
  
  sprintf(filename,"isodata");
  
  for (id=0;id<=MAX_ID;id++)
    {
      for (k=0;k<FS;k++)
	{
	  
	  sprintf(filename2,"%s.%d.%d.dat",filename,k,id); 
	  
	  /* Open tfile */
	  if (( fp = fopen(filename2, "r")) == NULL)
	    {
	      fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename2);
	      exit(0);
	    }

	  fscanf(fp,"%s",&junk);
	
	  for (i=0; i<=NBBINS; i++)
	    {	      
		  fscanf(fp,"%lf %lf \n",&junk2, &((*bbin_ptr).intensity[i][k][id]) );   		
		  
	    }
	  
	  fclose(fp); 

	}
    }
 

  
  return;
  
}
   
/****************** write cells*********************************/

void write_star_cells(s,c, ncells1,ncells2)

     int s;
     struct cell_type c[];
     int ncells1,ncells2;
     
{
  
  int i;
  double tmass=0.;
  char filename[20]; 
  FILE *fp;
  

  sprintf(filename,"incells.star.%d.dat",s); 
  
  /* Open tfile */
  if (( fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN cells.stardat FILE. \n");
      exit(0);
    }
  
  fprintf(fp, "@    x             y             z       nabs nscat  abslum      temp    dens       mass       cell_tau (normalised)\n");
  
  for (i=ncells1;i<=ncells2;i++)
    {
      fprintf(fp, " % .5e % .5e % .5e %4ld %4ld % .4e %6.1f % .4e % .4e % .4e %4ld\n", 
	      
	      c[i].x,
	      c[i].y, 
	      c[i].z,     
	      c[i].nabs,
	      c[i].nscat,
	      c[i].abslum,
	      c[i].temp,
	      c[i].dens,
	      c[i].mass,
	      c[i].tau,c[i].next_n); 
      
      tmass=tmass+c[i].mass ;
            
    }
  fclose(fp);
 
  
  printf ("total star.grid mass : %.4e Msun\n",tmass/SOLAR_MASS);
  
  return;
}



/***** function read_ptable *****************************************/

void read_ptable(reemit_ptr)

    struct reemit_type  *reemit_ptr;

{
    int i,k,flag;
    FILE *fp;

   /* Open  file */
   if (( fp = fopen("ptable.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ptable.dat FILE. \n");
         exit(0);
         }

   for (i = 1;i<TEMPDIM; i++)
     {
       flag=0;

       for (k=0; k<NFREQ ; k++)
	 {   
	   fscanf(fp,"%lf", & (*reemit_ptr).prob[i][k]); 
	   if (flag==0 && (*reemit_ptr).prob[i][k]==1.) {(*reemit_ptr).last[i]=k; flag=1;}
	 }

       if (flag==0) (*reemit_ptr).last[i]=NFREQ-1;
   }   
 
   fclose(fp);
   return;
}

/***** function read_rtable *****************************************/

void read_rtable(reemit_ptr)
     
     struct reemit_type  *reemit_ptr;

{
  int i,k,flag;
  FILE *fp;
  
  /* Open  file */
  if (( fp = fopen("rtable.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN rtable.dat FILE. \n");
      exit(0);
    }
  
  for (i = 1;i<TEMPDIM; i++)
    fscanf(fp,"%lf %d", & (*reemit_ptr).rfreq[i], & (*reemit_ptr).fpoint[i]);
   
  fclose(fp);
  
  return;
}



/***** function read_ptable *****************************************/

void read_ptable2(reemit2_ptr)

    struct reemit_type  *reemit2_ptr;

{
    int i,k,flag;
   FILE *fp;

 /* Open  file */
   if (( fp = fopen("ptable2.dat", "r")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN ptable.dat FILE. \n");
         exit(0);
         }

   for (i = 1;i<TEMPDIM; i++)
     {
       flag=0;

       for (k=0; k<NFREQ ; k++)
         {
           fscanf(fp,"%lf", & (*reemit2_ptr).prob[i][k]);
           if (flag==0 && (*reemit2_ptr).prob[i][k]==1.) {(*reemit2_ptr).last[i]=k; flag=1;}
         }

       if (flag==0) (*reemit2_ptr).last[i]=NFREQ-1;
   }

   fclose(fp);

   return;
}
/*********** function write direct photons ****************************************/

/* creates a file that has info about photons that directly escape from the system*/
/* for checking purporses */

void write_direct_photons(p,i,last_tau)

 struct photon_packet_type p[];
 long int  i;
 double last_tau;
 
{ 
FILE *fp;

   /* Open tfile */
   if (( fp = fopen("outdirect.dat", "a+")) == NULL)
         {
         fprintf(stderr, "FAILED TO OPEN outdirect.dat FILE. \n");
         exit(0);
         }

   fprintf(fp, "%ld % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e %.3f %5ld %5ld %5ld %2d\n",i, 
              p[i].x,
              p[i].y,
              p[i].z,
              p[i].kx,
              p[i].ky,
              p[i].kz,
              p[i].f,
              p[i].lum,
              last_tau,
	      p[i].starlum,
	      p[i].abs,
              p[i].scat,
              p[i].lcell,
              p[i].id
               );               
     
   fclose(fp);

}

/**************** write_restart_info ********************************/

void write_restart_info(s,output_counter,k_fpoints,photon_i,xnum,interactions,time0,rtime0)
     int s;
     double output_counter;
     long k_fpoints;
     double photon_i;
     long xnum;
     double interactions;
     double time0;
     double rtime0;
{
  
  char filename[30];
  FILE *fp;
  
  sprintf(filename,"restart.dat",s);
  
  if (( fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN restart.dat FILE. \n");
      exit(0);
    }
  
  fprintf(fp,"@#source / output_counter/ k_fpoints/ photon_i/xnum/ interactions / cpu time (min) \n");
  
  fprintf(fp,"%d %.0f %ld %.0f %ld %10.0f  %.2f\n",
	  s,output_counter,k_fpoints,photon_i,xnum,interactions,rtime0+(time(NULL)-time0)/60);
  
  fclose(fp);
  
  return;
}

/**************** read_restart_info ********************************/

void read_restart_info(s_start_ptr,output_counter_start_ptr,k_fpoints_start_ptr,photon_i_start_ptr,xnum_ptr,interactions_ptr,rtime0_ptr)
     int *s_start_ptr;
     double *output_counter_start_ptr;
     long *k_fpoints_start_ptr;
     double *photon_i_start_ptr;
     long *xnum_ptr;
     double *interactions_ptr;
     double *rtime0_ptr;

{
  char buff[MAX_LINE];
  int i;
  FILE *fp;
  
  if (( fp = fopen("restart.dat", "r")) == NULL)
     {
      fprintf(stderr, "FAILED TO OPEN restart.dat FILE. \n");
      exit(0);
    }
  
  for (i = 0; (i <5) & getdata(fp, buff); i++)
    { 
      sscanf(buff,"%d %lf %ld %lf %ld %lf %lf",
	     &(*s_start_ptr),&(*output_counter_start_ptr),
	     &(*k_fpoints_start_ptr),
	     &(*photon_i_start_ptr),&(*xnum_ptr),&(*interactions_ptr),&(*rtime0_ptr));
    }
  
  
  fclose(fp);
   
  return;
  
}


 
/**************** read interaction table ********************************/

void read_interactions(intercount_ptr,s,packets)

     double *intercount_ptr;
     int s;
     double packets;
{
  int t;char filename[30]; double value; int junk; char buff[MAX_LINE];
  FILE *fp;
  
  sprintf(filename,"intercount.%d.dat",s);

  if (( fp = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN %s FILE. \n",filename);
      exit(0);
    }

   for (t = 0; (t <60) & getdata(fp, buff); t++)
     { 
       sscanf(buff," %d %lf ",&junk,&value); 
       *(intercount_ptr+t)=value*packets/100.;
     }

  fclose(fp);
  
  return;
}
/****************** read outcells *********************************/

void read_active_cells_info(c, ncells_ptr)

 struct cell_type c[];
 int *ncells_ptr;

{
  char buff[MAX_LINE];
  int i,id;
  FILE *fp;
  
  
  
  if (( fp = fopen("outcells.sph.dat", "r")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN outcells.dat FILE. \n");
      exit(0);
    }
  
  for (i = 0; (i <MAX_CELLS) & getdata(fp, buff); i++)
    { 
     
      sscanf(buff,"%lf %lf %lf %ld %ld %lf %lf %lf %lf %d  %lf  %lf",
	     &(c[i].x),
	     &(c[i].y),
	     &(c[i].z),                                                  
	     &(c[i].nabs),
	     &(c[i].nscat),
	     &(c[i].abslum),
	     &(c[i].temp),
	     &(c[i].dens),
	     &(c[i].mass),
	     &id,
	     &(c[i].nsph),
	     &(c[i].size));         
      *ncells_ptr=(*ncells_ptr)+1;
    }
   
       

  fclose(fp);
 
  return;

}




/****************** read SPH DRAGON FILE  *********************************/

void rw_dragon_file(infilename,star,mass_ptr,nstar_ptr,dragon_time_ptr,dragon_step_ptr)

     char infilename[30];
     struct star_type star[];
     double *mass_ptr;
     int *nstar_ptr;
     double *dragon_time_ptr;
     long  *dragon_step_ptr;
 
{ int s=0;
  long c,nsph,nsph_tot;
  long i_dragon_info[20];
  double r_dragon_info[50];
  char outfilename[30];
  struct particle_type ptype;
  double tmass=0.;
  long nmass1,nmass2,nmass3;

  FILE *fp;

  ptypeptr pa=(ptypeptr)malloc(sizeof(ptype)*6000000);
  
  printf("\n--------------  Reading SPH DRAGON file --------------------------   \n");

  /*----*/
  if ((fp=fopen(infilename,"r")) == NULL)
    {
      fprintf(stderr, "initialize_from_out: FAILED TO OPEN %s file. \n",infilename);
      exit(0);
    }
  
  for (c = 0; c <20 ; c++) fscanf(fp,"%ld\n",&i_dragon_info[c]);

  *dragon_step_ptr=(long)(i_dragon_info[1]);  

  fscanf(fp,"%lf\n",dragon_time_ptr);  

  for (c = 0; c <49 ; c++) fscanf(fp,"%lf\n",&r_dragon_info[c]);
 
  nsph_tot=(long)(i_dragon_info[0]);
 

  printf("Info read ...\n");
  for (c = 0; c <=nsph_tot; c++) 
  { 
  fscanf(fp, "%lf %lf %lf", &(pa[c].x), &(pa[c].y), &(pa[c].z));
  }
  
    printf("Positions done... \n");
    
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf %lf %lf", &(pa[c].vx), &(pa[c].vy), &(pa[c].vz));   
 printf("Velocities done... \n");
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].temp));   
 printf("Temps done... \n");
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].h));   
 printf("h s done... \n");
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].rho));
 printf("Densities done... \n");
  for (c = 0; c <nsph_tot; c++) {fscanf(fp, "%lf", &(pa[c].mass)); }
 printf("Masses done... \n");
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%d", &(pa[c].iwas));
   printf("Itypes done... \n");
  fclose(fp);
  /*----*/

#ifdef RSMITH 

  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas==0) pa[c].iwas=1;

#endif 

  nsph=0; 

  for (c = 0; c <nsph_tot; c++) {if (pa[c].iwas>0) nsph=nsph+1;}
  /* calculate total mass */
  
  tmass=0.; 

  for (c = 0; c <nsph_tot; c++) {if (pa[c].iwas>0) tmass=tmass+ pa[c].mass; }
 
  *mass_ptr=tmass;
  

 printf("\nTIME........................... %.3e (Myr) \n",*dragon_time_ptr);
printf("\nTotal number of SPH particles.. %ld\n",i_dragon_info[0]);
 printf("\nSPH step ...................... %ld\n",*dragon_step_ptr);
 printf("Active SPH Particles........... %.0f\n",1.*nsph);  
 printf("sink radius.................... %.2e (AU)\n",r_dragon_info[37]*2e5*AU); 
 printf("Total active mass.............. %.4e (Msun)\n",tmass);


 sprintf(outfilename,"phaethon.in.%s.dat",infilename);
 
 
 /*----*/
 if (( fp = fopen(outfilename, "w")) == NULL)  
    {
      fprintf(stderr, "i. FAILED TO OPEN %s FILE. \n",outfilename);  
      exit(0);
    }

  /* WRITE IN BH FORMAT */

  fprintf(fp,"%.0f\n3 \n 5\n",1.*nsph);  
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .6e\n", pa[c].mass*SOLAR_MASS);
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .10e % .10e % .10e\n",pa[c].x*206265*AU,pa[c].y*206265*AU,pa[c].z*206265*AU);
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .10e % .10e % .10e\n",pa[c].vx,pa[c].vy,pa[c].vz);

  fclose(fp);



  /*----*/
  /* for (c = 0; c <nsph_tot; c++) if (pa[c].iwas<0) printf("%0.f\n",pa[c].iwas); */

 sprintf(outfilename,"phaethon.in.stars.%s.dat",infilename);

  if (( fp = fopen(outfilename, "w")) == NULL)
    {
      fprintf(stderr, "ii. FAILED TO OPEN %sFILE. \n",outfilename);
      exit(0);
    }
 
  /* WRITE STAR INFO */
  
  fprintf(fp,"%.0f\n3 \n 5\n",1.*nsph_tot);  
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .6e % .6e\n",*dragon_time_ptr,pa[c].mass*SOLAR_MASS);} 
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].x* 2e5*AU,pa[c].y*2e5*AU,pa[c].z*2e5*AU);}
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].vx,pa[c].vy,pa[c].vz);}
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .6e\n", pa[c].rho*SOLAR_MASS/pow(2e5*AU,3));}  

  fclose(fp);
  /*----*/


  for (c = 0; c <nsph_tot; c++) 
    {    
      if (pa[c].iwas<0)
	{ 
	  star[s].type=0; 
	  *nstar_ptr=(*nstar_ptr)+1;
	  /* printf("nstar=%d\n",*nstar_ptr);*/
	  star[s].mass=pa[c].mass*SOLAR_MASS;
	  star[s].x=pa[c].x* 2e5*AU;
	  star[s].y=pa[c].y* 2e5*AU;
	  star[s].z=pa[c].z* 2e5*AU;
	  printf("\n---------- STAR %d  -----------\n",s);
	  printf("Star mass: %.4e (Msun)\n",star[s].mass/SOLAR_MASS);
	  printf("Star x   :% .4e (AU)\n",star[s].x/AU   );
	  printf("Star y   :% .4e (AU)\n",star[s].y/AU   );
	  printf("Star z   :% .4e (AU)\n",star[s].z/AU   );
	  printf("Star x   :% .15e (cm)\n",star[s].x  );
	  printf("Star y   :% .15e (cm)\n",star[s].y  );
	  printf("Star z   :% .15e (cm)\n",star[s].z  );
	  s=s+1; if (s>10) printf("Error: too many stars\n");
	}

    }

  if (ISRF==1) 
    {

      printf("check this out inout.c -- F00\n");
	
      exit(0);
      *nstar_ptr=(*nstar_ptr)+1;

      star[(*nstar_ptr)].type=1;

      printf("\nIncluding ISRF................. \n");
    }


  if (MANUAL_STAR==1)
    {
      printf("check this out inout.c -- F01\n");
      exit(0);
      *nstar_ptr=(*nstar_ptr)+1;
      
      star[(*nstar_ptr)].type=MANUAL_STAR_TYPE;
      
      printf("\nIncluding MANUAL_STAR........... \n");
    }
 
  printf("NUMBER OF RADIATION SOURCES.... %d \n",(*nstar_ptr)+1);
  free(pa);
  pa=NULL;
  
  printf("\n--------------   END  reading SPH DRAGON file ------------------------ --   \n");

 
  return;
}
   


/****************** read SPH DRAGON FILE  *********************************/

void rw_spyros_file(infilename,star,mass_ptr,nstar_ptr,dragon_time_ptr,dragon_step_ptr)

     char infilename[30];
     struct star_type star[];
     double *mass_ptr;
     int *nstar_ptr;
     double *dragon_time_ptr;
     long  *dragon_step_ptr;
 
{ int s=0;
  int c,nsph,nsph_tot;
  int i_dragon_info[20];
  double r_dragon_info[50];
  char outfilename[30];
  struct particle_type ptype;
  double tmass=0.;
  long nmass1,nmass2,nmass3;

  FILE *fp;

  ptypeptr pa=(ptypeptr)malloc(sizeof(ptype)*325000);
  
  printf("\n--------------  Reading SPH SPYROS file --------------------------   \n");

  /*----*/
  if ((fp=fopen(infilename,"r")) == NULL)
    {
      fprintf(stderr, "initialize_from_out: FAILED TO OPEN %s file. \n",infilename);
      exit(0);
    }
  
  
 
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf %lf %lf", &(pa[c].x), &(pa[c].y), &(pa[c].z));
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf %lf %lf", &(pa[c].vx), &(pa[c].vy), &(pa[c].vz));   
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].temp));   
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].h));   
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].rho));
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%lf", &(pa[c].mass));
  for (c = 0; c <nsph_tot; c++) fscanf(fp, "%d", &(pa[c].iwas));
 
  fclose(fp);
  /*----*/


  nsph=0; 

  for (c = 0; c <nsph_tot; c++) {if (pa[c].iwas>0) nsph=nsph+1;}
  /* calculate total mass */
  
  tmass=0.; 

  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) tmass=tmass+ pa[c].mass;  
 
  *mass_ptr=tmass;
  

 printf("\nTIME........................... %.3e (Myr) \n",*dragon_time_ptr);
 printf("\nTotal number of SPH particles.. %d\n",i_dragon_info[0]);
 printf("\nSPH step ...................... %ld\n",*dragon_step_ptr);
 printf("Active SPH Particles........... %.0f\n",1.*nsph);  
 printf("sink radius.................... %.2e (AU)\n",r_dragon_info[37]*2e5*AU); 
 printf("Total active mass.............. %.4e (Msun)\n",tmass);


 sprintf(outfilename,"phaethon.in.%s.dat",infilename);
 
 
 /*----*/
 if (( fp = fopen(outfilename, "w")) == NULL)  
    {
      fprintf(stderr, "FAILED TO OPEN %s FILE. \n",outfilename);  
      exit(0);
    }

  /* WRITE IN BH FORMAT */

  fprintf(fp,"%.0f\n3 \n 5\n",1.*nsph);  
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .6e\n", pa[c].mass*SOLAR_MASS);
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].x*2e5*AU,pa[c].y*2e5*AU,pa[c].z*2e5*AU);
  for (c = 0; c <nsph_tot; c++) if (pa[c].iwas>0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].vx,pa[c].vy,pa[c].vz);

  fclose(fp);



  /*----*/
  /* for (c = 0; c <nsph_tot; c++) if (pa[c].iwas<0) printf("%0.f\n",pa[c].iwas); */

 sprintf(outfilename,"phaethon.in.stars.%s.dat",infilename);

  if (( fp = fopen(outfilename, "w")) == NULL)
    {
      fprintf(stderr, "FAILED TO OPEN %sFILE. \n",outfilename);
      exit(0);
    }
 
  /* WRITE STAR INFO */
  
  fprintf(fp,"%.0f\n3 \n 5\n",1.*nsph_tot);  
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .6e % .6e\n",*dragon_time_ptr,pa[c].mass*SOLAR_MASS);} 
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].x* 2e5*AU,pa[c].y*2e5*AU,pa[c].z*2e5*AU);}
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .4e % .4e % .4e\n",pa[c].vx,pa[c].vy,pa[c].vz);}
  for (c = 0; c <nsph_tot; c++) 
    {if (pa[c].iwas<0) fprintf(fp,"% .6e\n", pa[c].rho*SOLAR_MASS/pow(2e5*AU,3));}  

  fclose(fp);
  /*----*/


  for (c = 0; c <nsph_tot; c++) 
    {    
      if (pa[c].iwas<0)
	{ 
	  star[s].type=0; 
	  *nstar_ptr=(*nstar_ptr)+1;
	  /* printf("nstar=%d\n",*nstar_ptr);*/
	  star[s].mass=pa[c].mass*SOLAR_MASS;
	  star[s].x=pa[c].x* 2e5*AU;
	  star[s].y=pa[c].y* 2e5*AU;
	  star[s].z=pa[c].z* 2e5*AU;
	  printf("\n---------- STAR %d  -----------\n",s);
	  printf("Star mass: %.4e (Msun)\n",star[s].mass/SOLAR_MASS);
	  printf("Star x   :% .4e (AU)\n",star[s].x/AU   );
	  printf("Star y   :% .4e (AU)\n",star[s].y/AU   );
	  printf("Star z   :% .4e (AU)\n",star[s].z/AU   );
	  printf("Star x   :% .15e (cm)\n",star[s].x  );
	  printf("Star y   :% .15e (cm)\n",star[s].y  );
	  printf("Star z   :% .15e (cm)\n",star[s].z  );
	  s=s+1; if (s>10) printf("Error: too many stars\n");
	}

    }

  if (ISRF==1) 
    {
      *nstar_ptr=(*nstar_ptr)+1;

      star[(*nstar_ptr)].type=1;

      printf("\nIncluding ISRF................. \n");
    }



  printf("NUMBER OF RADIATION SOURCES.... %d \n",(*nstar_ptr)+1);
  free(pa);
  pa=NULL;
  
  return;
}
   
