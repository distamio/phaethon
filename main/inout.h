void read_treebo(ptypeptr,long int  *);
void write_treebo(ptypeptr,long int);
void write_photons(phtypeptr);
void read_photons(phtypeptr);
void write_cells(ctypeptr, int);
void write_star_cells(int,ctypeptr,int,int);
void read_cells(ctypeptr, long *);
void read_cells_xyz(ctypeptr, long *);
void write_grid(struct grid_type);
void read_grid(struct grid_type *);
void read_opacity(struct opacity_type *);
void write_opacity(struct opacity_type);
void write_mean_opacity(struct mean_op_type);
void read_mean_opacity(struct mean_op_type *);
void record_photons(phtypeptr);
void record_cells(ctypeptr, int);
void read_outcells(ctypeptr, long *);
void read_outphotons(phtypeptr);
void record_sed_info(phtypeptr);
void read_sed_info(phtypeptr);
void write_ftable(struct reemit_type);
void read_ftable(struct reemit_type *);
void read_ptable(struct reemit_type *);
void read_rtable(struct reemit_type *);
void read_ISRF(struct btype *);
void record_photon_basic(phtypeptr,double);
void read_outphotons_basic(phtypeptr);
void read_isodata(struct bbin_type *, int, int); 
void read_bbins(struct bbin_type *);
void record_active_cells(ctypeptr, int);
void read_active_cells_info(ctypeptr, int *);
void   write_direct_photons(phtypeptr,long,double);
void write_info_file(phtypeptr);
void record_potential_cell_info(ctypeptr, int,double);
void read_active_cells(ctypeptr, int);
void record_active_cell_info(ctypeptr, int,double);
void rw_dragon_file();
void rw_GCODE_file();
void record_active_cell_info_POSTRT(ctypeptr, int, double);
void read_star_sed(startypeptr,int);
 

#ifdef DUST_MIX
void read_opacity2(struct opacity_type *);
void read_mean_opacity2(struct mean_op_type *);
void read_ftable2(struct reemit_type *);
void read_ptable2(struct reemit_type *);
#endif


#ifdef DO_ISOMAP
void read_isogrid(struct isogrid_type *, int, int);
#endif


void write_restart_info(int,double, long, double,long,double,double,double);
void read_restart_info(int *,double *, long *, double *,long *,double*,double *);
void read_interactions(double *, int, double);
