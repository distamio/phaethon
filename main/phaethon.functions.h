#include <string.h>
void   advance_photon    (phtypeptr,int,ctypeptr,int,startypeptr,unsigned short);
void   advance_photon_prt(phtypeptr,int,ctypeptr,int, startypeptr,int);
void   advance_photon_free(phtypeptr,int,ctypeptr,int,double,startypeptr,unsigned short);
void   advance_external_photon(phtypeptr,long,ctypeptr,double,double);

int    find_photon_cell(phtypeptr,long,ctypeptr,startypeptr,int);
int    find_photon_cell_advance(phtypeptr,long,ctypeptr,startypeptr,int);
int    find_prt_photon_cell(phtypeptr,long,ctypeptr,startypeptr,int);
int    find_test_photon_cell(phtypeptr,long,ctypeptr,startypeptr,int);

double find_mfp_free(ctypeptr,int);

double find_stepsize(phtypeptr,long,ctypeptr);
int    find_cell_from_position(double,ctypeptr);

void   find_photon_angles(phtypeptr,int,ctypeptr);
void   find_photon_dust_kind(phtypeptr,int);
void   find_photon_origin(phtypeptr,int,startypeptr,int);  
void   make_prt_photon(phtypeptr,ctypeptr,long ,long,double,double, ctypeptr, long,long,startypeptr,int,int);  

void   make_test_photon(phtypeptr,long,double,double,double,double,double);
void   advance_test_photon(phtypeptr,long,ctypeptr,int,double,double,startypeptr, int);
/* void   make_cells(ctypeptr,int *); */


