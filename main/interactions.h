

double ree(double, double, double, struct mean_op_type *);
double solve_ree(double, double, double,int,struct mean_op_type *); 
double   absorb_photon(double, long, ctypeptr,struct mean_op_type *);
void   reemit_photon(phtypeptr, long,struct reemit_type *);
void   find_new_f(phtypeptr,long,struct reemit_type *);
void   fix_temp(ctypeptr);
int invree(double, double, double, struct mean_op_type *);
void   scatter_photon(phtypeptr,long);
#ifndef DUST_MIX
void  find_cells_next_n(ctypeptr, long, double, struct mean_op_type *);
#endif

#ifdef DUST_MIX

double ree2(double, double, double, struct mean_op_type *);
double solve_ree2(double, double, double,int,struct mean_op_type *);
double   absorb_photon2(double, long, ctypeptr,struct mean_op_type *);
void   reemit_photon2(phtypeptr, long,struct reemit_type *);
void   find_new_f2(phtypeptr,long,struct reemit_type *);
int invree2(double, double, double, struct mean_op_type *);
void  find_cells_next_n(ctypeptr, long, double, struct mean_op_type *,struct mean_op_type *);
#endif


