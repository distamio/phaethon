double k_abs(double,struct opacity_type);
double sigma_scat(double,struct opacity_type);
double albedo(double,struct opacity_type);
double mean_opacity(double,struct mean_op_type *);
int binary_search_opacity(double, int,int,struct opacity_type);
void  create_opacity_table (struct opacity_table *, struct reemit_type  *, struct opacity_type);
void   find_photon_fpoint(phtypeptr,int,struct reemit_type *);
