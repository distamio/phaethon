#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "constants.h"
#include "phaethon.params.h"
#include "params.grid.h"
#include  "params.h"

void read_params ()
{
char line[100];
FILE *fp;

if ((fp=fopen("params.dat","r"))==NULL)
     {
      fprintf(stderr, "FAILED TO OPEN params.grid.dat FILE. \n");
      exit(0);
      }

// CONTROL OPTIONS 
fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%lf\n",&OUTPUT_FACTOR); fgets(line, 100, fp);
fscanf(fp,"%d\n",&RESTART); fgets(line, 100, fp);
fscanf(fp,"%d\n",&MORE); fgets(line, 100, fp);

// AMBIENT RADIATION FIELD (ISRF)

fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%d\n",&ISRF); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&ISRF_NPACKETS); fgets(line, 100, fp);
fscanf(fp,"%d\n",&R_ISRF_R_MAX); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&R_ISRF); fgets(line, 100, fp);
R_ISRF=R_ISRF*LENGTH_UNIT;
if (R_ISRF_R_MAX==1) R_ISRF=R_MAX;
fscanf(fp,"%lf\n",&ISRF_MULTIPLY); fgets(line, 100, fp);
fscanf(fp,"%d\n",&ISRF_IGNORE_SCATTER); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STAR_ISRF_TEMP); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&DILUTION); fgets(line, 100, fp);

// STELLAR RADIATION FIELD (STAR)  

fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%d\n",&MANUAL_STAR); fgets(line, 100, fp);
fscanf(fp,"%d\n",&MANUAL_STAR_TYPE); fgets(line, 100, fp);
fscanf(fp,"%d\n",&STAR_SED_FILE); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STAR_NPACKETS); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STAR_RADIUS); fgets(line, 100, fp);
STAR_RADIUS=STAR_RADIUS*SOLAR_RADIUS;
fscanf(fp,"%lf\n",&STAR_X); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STAR_Y); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STAR_Z); fgets(line, 100, fp);
STAR_X=STAR_X*LENGTH_UNIT;
STAR_Y=STAR_Y*LENGTH_UNIT;
STAR_Z=STAR_Z*LENGTH_UNIT;
fscanf(fp,"%lf\n",&STAR_MASS); fgets(line, 100, fp);
STAR_MASS=STAR_MASS*SOLAR_MASS;
fscanf(fp,"%lf\n",&STAR_TEMP); fgets(line, 100, fp);

//SED & ISOPHOTAL MAPS PARAMETERs
fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%lf\n",&DISTANCE); fgets(line, 100, fp);
DISTANCE=DISTANCE*PC;
fscanf(fp,"%lf\n",&DTHETA); fgets(line, 100, fp);
DTHETA=DTHETA*M_PI;
fscanf(fp,"%lf\n",&SL1); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&SL2); fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNFBINS); fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNBBINS); fgets(line, 100, fp);
// Wavelengths (microns) to calculate isophotal maps (up to 10)
fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iFS); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F0,&dF0); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F1,&dF1); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F2,&dF2); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F3,&dF3); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F4,&dF4); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F5,&dF5); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F6,&dF6); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F7,&dF7); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F8,&dF8); fgets(line, 100, fp);
fscanf(fp,"%lf %lf \n",&F9,&dF9); fgets(line, 100, fp);

//  Number of pixels of intensity maps 
//fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNXPIXELS); fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNYPIXELS); fgets(line, 100, fp);
//Size of area of intensity map (pc or AU)
fgets(line, 100, fp);
fscanf(fp,"%lf\n",&XMIN); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&XMAX); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&YMIN); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&YMAX); fgets(line, 100, fp);
XMIN=XMIN*LENGTH_UNIT;
YMIN=YMIN*LENGTH_UNIT;
XMAX=XMAX*LENGTH_UNIT;
YMAX=YMAX*LENGTH_UNIT;
// Polar angles to compute intensity maps (up to 10)
fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNLBINS); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA0); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA1); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA2); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA3); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA4); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA5); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA6); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA7); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA8); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&THETA9); fgets(line, 100, fp);
//  Azimuthal angles to compute intensity maps (up to 10)
fgets(line, 100, fp);
//fscanf(fp,"%d\n",&iNMBINS); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI0); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI1); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI2); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI3); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI4); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI5); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI6); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI7); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI8); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&PHI9); fgets(line, 100, fp);
//  DUST PARAMETERS  
fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%lf\n",&DUST_TEMP); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&SCATTER_MULTIPLY); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&G_SCATTER); fgets(line, 100, fp);
// PROPAGE PHOTONS
fgets(line, 100, fp);fgets(line, 100, fp);fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STEP_MFP); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STEP_DENS); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STEP_R); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&STEP_CELL); fgets(line, 100, fp);
fscanf(fp,"%lf\n",&INNER_STEP_SIZE); fgets(line, 100, fp);
fclose(fp);

return;
 }
