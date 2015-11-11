void init_isrf_photons_3D(long,phtypeptr,startypeptr,int);
void init_isrf_photons_1D(long,phtypeptr,startypeptr,int);
void init_star_photons( long,phtypeptr, startypeptr, int);
void find_xtable(struct xtable_type *,startypeptr,int);
void assign_photon_freq_lum(long,int *,phtypeptr,struct xtable_type *,startypeptr, int,int *, long*);
void do_diagnostics(phtypeptr,startypeptr, int);
void find_xtable_isrf(struct xtable_type *,struct btype *,startypeptr, int);
void find_xtable_multiplicity(struct xtable_type *);
void init_photon_params_box(long,phtypeptr,double);
void assign_photon_opacity(phtypeptr,int,struct opacity_table *, int);
void assign_photon_opacity2(phtypeptr,int,struct opacity_table *, int);
void find_xtable_star_sed(struct xtable_type *,startypeptr,int);

