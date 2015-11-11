#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "params.grid.h"
#include "structs.h"
#include "opacity.h"
#ifdef DUST_MIX

void read_grid_params (opacity,opacity2)

 struct opacity_type opacity;
 struct opacity_type opacity2;

#else

void read_grid_params (opacity)

  struct opacity_type opacity;
#endif

{
char line[100];
FILE *fp;


if ((fp=fopen("params.grid.dat","r"))==NULL)
     {
      fprintf(stderr, "FAILED TO OPEN params.grid.dat FILE. \n");
      exit(0);
      }
printf("\n::Reading grid parameters (params.grid.dat)\n");

fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);

// Read grid parameters
fscanf(fp,"%d\n",&ENVELOPE); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&CDENS0); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R0); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R_MAX_CORE); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R_MAX); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R_ORIGIN); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&DENS_EXT_FACTOR); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R_DUST); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&MU); fgets(line, 100, fp);
fscanf(fp,"%d\n",&BECELLS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&MCCELLS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&NLCELLS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&NMCELLS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&LINEAR); fgets(line, 100, fp);
fscanf(fp,"%d\n",&ISOANGLES); fgets(line, 100, fp);
fscanf(fp,"%d\n",&FTAU_POINTS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&GEOMETRY); fgets(line, 100, fp);
fscanf(fp,"%d\n",&SIN_POWER); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&ASYMMETRY_FACTOR); fgets(line, 100, fp);

if (SIN_POWER!=1 && SIN_POWER!=2 && SIN_POWER!=4) 
{
 printf ("\n:: grid.params.c :: SIN_POWER:: Invalid value !!\n");
 exit(0);
}

if (ASYMMETRY_FACTOR>3 && ASYMMETRY_FACTOR<=1 ) 
{
 printf ("\n:: grid.params.c :: ASYMMETRY_FACTOR>3 or <=1 :: Invalid value !!\n");
 exit(0);
}

A=59*(ASYMMETRY_FACTOR-1.);

// Transform grid parameters
R0=R0*LENGTH_UNIT;
R_MAX_CORE=R_MAX_CORE*LENGTH_UNIT;
R_MAX=R_MAX*LENGTH_UNIT;
R_DUST=R_DUST*R_MAX_CORE;
R_ORIGIN=R_ORIGIN*LENGTH_UNIT;
if (MCCELLS==0) 
{
	printf("\n:: MCCELLS=0 :: Setting R_MAX=R_MAX_CORE");
	R_MAX=R_MAX_CORE;
}
if (R_MAX<R_MAX_CORE && MCCELLS<0)
{
	printf("\n:: ERROR: R_MAX<R_MAX_CORE");
	exit(0);
}


// Additional grid parameters
DENS0=CDENS0*MU*M_PROTON;
NRCELLS=BECELLS+MCCELLS;
NCELLS=NLCELLS*NMCELLS*NRCELLS;
K_VISUAL0=k_abs((c_CONST/0.55)*1e4,opacity)+sigma_scat((c_CONST/0.55)*1e4,opacity);
printf("\n:: VISUAL EXTINCTIONS :: K_VISUAL0= %.2e (cm2/g)\n", K_VISUAL0);

#ifdef DUST_MIX
K_VISUAL1=k_abs((c_CONST/0.55)*1e4,opacity2)+sigma_scat((c_CONST/0.55)*1e4,opacity2);
printf(":: VISUAL EXTINCTIONS :: K_VISUAL1= %.2e (cm2/g)\n\n", K_VISUAL1);
#endif

fscanf(fp,"%d\n",&FTAU_POINTS); fgets(line, 100, fp);

fgets(line, 70, fp); fgets(line, 70, fp); fgets(line, 70, fp);
fscanf(fp, "%d\n",&DISC);  fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_CDENS0); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_R); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_ALPHA); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_INNER_GAP); fgets(line, 100, fp);
fscanf(fp, "%d\n",&DISC_INNER_R_DUST); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_THETA_OPENING); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&DISC_AMBIENT_DENS); fgets(line, 100, fp);
fgets(line, 70, fp); fgets(line, 70, fp); fgets(line, 70, fp);
fscanf(fp, "%d\n",&JET);  fgets(line, 100, fp);
fscanf(fp, "%lf\n",&JET_DENS); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&JET_BETA); fgets(line, 100, fp);
fscanf(fp, "%lf\n",&JET_THETA); fgets(line, 100, fp);

if (JET==1) {
JET_THETA=JET_THETA*M_PI;
JET_ALPHA=(pow(R_MAX_CORE,1-JET_BETA)*cos(JET_THETA)*pow(sin(JET_THETA),-JET_BETA));
JET_DENS=JET_DENS*MU*M_PROTON;
}

if (DISC==1) 
{
DISC_R          = DISC_R*AU;
DISC_DENS0      = DISC_CDENS0*MU*M_PROTON;     
DISC_Z          = DISC_R*DISC_ALPHA;
DISC_INNER_GAP  = DISC_INNER_GAP*AU;
if (DISC_INNER_R_DUST==1) DISC_INNER_GAP=R_DUST;
if (DISC_INNER_GAP>R_DUST) R_DUST=DISC_INNER_GAP;
}
fclose(fp);

printf(":: DENS0=%.4e (g/cm3)\n",DENS0);
printf(":: tau_visual=K_VISUAL0*R_MAX_CORE*DENS0 : %.4e\n\n ", K_VISUAL0*R_MAX_CORE*DENS0);


#ifdef SPHERICAL

if (GEOMETRY>0) 
{
 printf ("\n:: grid.params.c :: SPHERICAL is defined and GEOMETRY is not 0 !!\n");
 exit(0);
}

if (DISC==1) 
{
 printf ("\n:: WARNING :: grid.params.c :: SPHERICAL is defined and DISC is 1 !!\n");
}
#endif

return;
 }
