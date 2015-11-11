
// CONTROL OPTIONS   
double OUTPUT_FACTOR;
int RESTART;
int MORE;
//AMBIENT RADIATION FIELD (ISRF) 
double ISRF_NPACKETS;
int ISRF;
int R_ISRF_R_MAX;
double R_ISRF;
double ISRF_MULTIPLY;
int ISRF_IGNORE_SCATTER;
double STAR_ISRF_TEMP;
double DILUTION;
// STELLAR RADIATION FIELD (STAR)  
double STAR_NPACKETS;
double STAR_RADIUS; 
double STAR_X;
double STAR_Y;
double STAR_Z;
double STAR_MASS;
double STAR_TEMP;
int MANUAL_STAR;
int MANUAL_STAR_TYPE;
int STAR_SED_FILE;
// SED & ISOPHOTAL MAPS PARAMETERs
double DISTANCE;
double DTHETA;
double SL1;
double SL2;
//int    iNFBINS;
//int    iNBBINS;
// Wavelengths (microns) to calculate isophotal maps (up to 10)
//int iFS;
double F0,F1,F2,F3,F4,F5,F6,F7,F8,F9;
double dF0,dF1,dF2,dF3,dF4,dF5,dF6,dF7,dF8,dF9;
//  Number of pixels of intensity maps
//int iNXPIXELS;
//int iNYPIXELS;
// Size of area of intensity map (pc or AU)
double XMIN,XMAX,YMIN,YMAX;
//  Polar angles to compute intensity maps (up to 10)
//int iNLBINS;
double THETA0,THETA1,THETA2,THETA3,THETA4,THETA5,THETA6,THETA7,THETA8,THETA9;
//Azimuthal angles to compute intensity maps (up to 10)
//int iNMBINS;
double PHI0,PHI1,PHI2,PHI3,PHI4,PHI5,PHI6,PHI7,PHI8,PHI9;
//  DUST PARAMETERS 
double DUST_TEMP;
double SCATTER_MULTIPLY;
double G_SCATTER;
//  PROPAGE PHOTONS  
double STEP_MFP;
double STEP_DENS;
double STEP_R;
double STEP_CELL;
double INNER_STEP_SIZE;
