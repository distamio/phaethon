
void create_fbins(struct fbin_type *,struct opacity_type);
void init_fbins(struct fbin_type *,struct opacity_type);
void place_into_fbins(phtypeptr,int,struct fbin_type *,double *,int,int,startypeptr,int);
void create_dbins(struct dbin_type *);
void init_dbins(struct dbin_type *);
void place_into_dbins(phtypeptr,int, struct dbin_type *, double *);
void find_sed(phtypeptr,struct fbin_type *,int,int,int,int);
void write_sed(struct fbin_type *,struct btype *,struct dbin_type,double,int *,int);
void write_sed_prt(struct fbin_type *,struct btype *,struct dbin_type,double,ctypeptr,ctypeptr,ctypeptr,startypeptr,int,int);
void create_bbins(struct bbin_type *,double);
void init_bbins(struct bbin_type *);
void place_into_bbins(phtypeptr,int, struct bbin_type *, double *,int*, int*, int*,double,startypeptr,int);
void write_isophotal_map2(struct bbin_type *,char[]);
int binary_search_bbin (double, int, int, struct bbin_type *);
void find_bbin_intensity(struct bbin_type *,double *,double *,startypeptr,int, int*);
int binary_search_fbin (double, int, int, struct fbin_type *);
void null_bbin_npackets (struct bbin_type *);

#ifdef DO_ISOMAP
int binary_search_isogrid_x (double, int, int, struct isogrid_type *,double);
int binary_search_isogrid_y (double, int, int, struct isogrid_type *,double);
void find_isogrid_intensity(struct isogrid_type *, double *,double *,struct dbin_type, double *,startypeptr,int,double,int*, int);
void write_isophotal_map_(struct isogrid_type);
void make_isogrid_spherical_case(struct bbin_type *,struct isogrid_type *,double *, int);
void make_isogrid(struct isogrid_type *);
void place_into_isogrid(phtypeptr,int,struct isogrid_type *,double, double,int*, int*, int*,int,int,startypeptr,int);
void calculate_isogrid_intensity(int, int, struct bbin_type *,struct isogrid_type *, double *, int);

void isomap_puma(struct isogrid_type *);
void isomap_idl(struct isogrid_type *);
void isomap_statistics_puma(struct isogrid_type *,int);
void isomap_statistics_idl(struct isogrid_type *,int);

void find_fbin_isogrid_3D(phtypeptr, int, struct dbin_type, struct isogrid_type *, struct fbin_type *, int *, int *, int* , double *,startypeptr,int);
void find_fbin_isogrid_2D_I(phtypeptr,int, struct dbin_type, struct isogrid_type *, struct fbin_type *,  int *, int*, int*, double *,startypeptr,int);
void find_fbin_isogrid_2D_II(phtypeptr,int, struct dbin_type, struct isogrid_type *, struct fbin_type *,  int *, int*, int *, double *,startypeptr,int); 

void null_isogrid_npackets(struct isogrid_type *);
			    
#endif

void calculate_domega(struct dbin_type *);
void find_obs_fbins(int *,struct fbin_type, double *, int *, int*);
