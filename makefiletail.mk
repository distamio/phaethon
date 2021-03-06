#makefiletail.mk
#
# Remove trailing whitespace from user options
# -----------------------------------------------------------------------------
# PHAETHON MAKEFILE 
# -----------------------------------------------------------------------------
CC               		 := $(strip $(CC))
OPENMP				 := $(strip $(OPENMP))
FASTMATH			 := $(strip $(FASTMATH))
OPTIMISE           	         := $(strip $(OPTIMISE))
PROFILE 			 := $(strip $(PROFILE))
MAIN_DIR             		 := $(strip $(MAIN_DIR))
GRID_DIR           		 := $(strip $(GRID_DIR))
OPACITY_DIR        		 := $(strip $(OPACITY_DIR))
OPACITY2_DIR                     := $(strip $(OPACITY2_DIR))
OPACITY_UTILS_DIR                := $(strip $(OPACITY_UTILS_DIR))
UTILS_DIR  			 := $(strip $(UTILS_DIR))
ISRF_DIR           		 := $(strip $(ISRF_DIR))
EXEC_DIR			 := $(strip $(EXEC_DIR))

# -----------------------------------------------------------------------------
# GRID 
# -----------------------------------------------------------------------------
SPHERICAL          	 	 := $(strip $(SPHERICAL))
DIM		            	 := $(strip $(DIM))

# -----------------------------------------------------------------------------
# SEDs  & INTENSITY MAPS
# -----------------------------------------------------------------------------
NFBINS				 := $(strip $(NFBINS))
FS				 := $(strip $(FS))
NLBINS				 := $(strip $(NLBINS))
NMBINS				 := $(strip $(NMBINS))

NBBINS				 := $(strip $(NBBINS))
NXPIXELS 			 := $(strip $(NXPIXELS))
NYPIXELS 			 := $(strip $(NYPIXELS))

ISODIM		    		 := $(strip $(ISODIM))	
EXTRA_SYMMETRY 	 		:= $(strip $(EXTRA_SYMMETRY))
SPHERICAL_ISOMAP   		:= $(strip $(SPHERICAL_ISOMAP))
SPHERICAL_SED    		:= $(strip $(SPHERICAL_SED))

# -----------------------------------------------------------------------------
# DUST PROPERTIES
# -----------------------------------------------------------------------------
DUST_MIX                         := $(strip $(DUST_MIX))
DUST_KIND                        := $(strip $(DUST_KIND))

# -----------------------------------------------------------------------------
# CONTROL OPTIONS 
# -----------------------------------------------------------------------------
NO_INTERACTION     	 := $(strip $(NO_INTERACTION))
MULTIPLE_PHOTONS   	 := $(strip $(MULTIPLE_PHOTONS))
TREE_INFO_OUT      	 := $(strip $(TREE_INFO_OUT))
INFO_OUTPUT        	 := $(strip $(INFO_OUTPUT))
INC_SHADOW_PHOTONS 	 := $(strip $(INC_SHADOW_PHOTONS))
DO_ISOMAP	   	 := $(strip $(DO_ISOMAP))
LENGTH_UNIT		 :=$(strip $(LENGTH_UNIT))
SPH_SIMULATION           := $(strip $(SPH_SIMULATION))

# Precision.  Possible values are SINGLEPREC, MIXEDPREC, and DOUBLEPREC (for tree)
PREC = SINGLEPREC

MAIN_HEADERS=$(MAIN_DIR)
GRID_HEADERS=$(GRID_DIR)
OPACITY_HEADERS+=$(OPACITY_DIR)

VPATH = $(MAIN_DIR):$(GRID_DIR):$(OPACITY_DIR):$(OPACITY_UTILS_DIR):$(ISRF_DIR):$(UTILS_DIR)

# Object files always included in compilation list
# ----------------------------------------------------------------------------
OBJ += dens.o advance.o interactions.o mgrid.o photons.o bb.o diagnostics.o grid.params.o opacity.o spec.o find_cell.o inout.o params.o star.o
OBJ_OPA=  opacity.o inout.o bb.o 

# Compiler flags
# ----------------------------------------------------------------------------
ifeq ($(CC),gcc)
OPT += -I $(MAIN_HEADERS)  -I $(GRID_HEADERS) -I $(OPACITY_HEADERS)  
endif
ifeq ($(CC),cc)
OPT += -I $(MAIN_HEADERS)  -I $(GRID_HEADERS) -I $(OPACITY_HEADERS)  
endif

ifeq ($(FASTMATH),1)
OPT += -ffast-math  -msse2 -mmmx
endif

ifeq ($(OPENMP),1)
#if is not working try fopenmp or mp
OPT += -openmp
endif

# Optimisation level
# ----------------------------------------------------------------------------
ifeq ($(OPTIMISE),1)
OPT += -O1
endif
ifeq ($(OPTIMISE),2)
OPT += -O2
endif
ifeq ($(OPTIMISE),3)
OPT += -O3  
endif
ifeq ($(OPTIMISE),4)
OPT += -O4
endif
ifeq ($(OPTIMISE),5)
OPT += -O5
endif

# Profiling options
# ----------------------------------------------------------------------------
ifeq ($(PROFILE),1)
OPT += -pg
endif
ifeq ($(PROFILE),2)
OPT += -pg -g
endif


# Grid flags
# ----------------------------------------------------------------------------
ifeq ($(SPHERICAL),1)
GFLAGS += -DSPHERICAL
endif

ifeq ($(DIM),4)
GFLAGS += -DDIM=4
endif
ifeq ($(DIM),2)
GFLAGS += -DDIM=2
endif
ifeq ($(DIM),3)
GFLAGS += -DDIM=3
endif

# SED &  isomaps flangs
# ----------------------------------------------------------------------------

SFLAGS +=-DNFBINS=$(NFBINS)
SFLAGS +=-DFS=$(FS)
SFLAGS +=-DNLBINS=$(NLBINS)
SFLAGS +=-DNMBINS=$(NMBINS)

ifeq ($(DO_ISOMAP),1)
SFLAGS += -DDO_ISOMAP=1
endif

SFLAGS +=-DNBBINS=$(NBBINS)
SFLAGS +=-DNXPIXELS=$(NXPIXELS)
SFLAGS +=-DNYPIXELS=$(NYPIXELS)

SFLAGS += -DISODIM=$(ISODIM)

SFLAGS +=-DEXTRA_SYMMETRY=$(EXTRA_SYMMETRY)

SFLAGS +=-DSPHERICAL_ISOMAP=$(SPHERICAL_ISOMAP)


ifeq ($(SPHERICAL_SED),1)
SFLAGS +=-DSPHERICAL_SED
endif

# SED &  isomaps flangs
# ----------------------------------------------------------------------------

ifeq ($(DUST_MIX),1)
OFLAGS += -DDUST_MIX
OBJ += dustkind.o 

ifeq ($(DUST_KIND),1)
OFLAGS += -DDISC_IN_CLOUD
endif

endif

# Other options
# ----------------------------------------------------------------------------
ifeq ($(NO_INTERACTION),1)
OFLAGS += -DNO_INTERACTION=1
endif

OFLAGS += -DMULTIPLE_PHOTONS=$(MULTIPLE_PHOTONS)

ifeq ($(TREE_INFO_OUT),1)
OFLAGS += -DTREE_INFO_OUT=1
endif

ifeq ($(INFO_OUTPUT),1)
OFLAGS += -DINFO_OUTPUT=1
endif

ifeq ($(INC_SHADOW_PHOTONS),1)
OFLAGS += -DINC_SHADOW_PHOTONS=1
endif

ifeq ($(LENGTH_UNIT),PC)
OFLAGS +=-DLENGTH_UNIT=3.09e18
endif

ifeq ($(LENGTH_UNIT),AU)
OFLAGS +=-DLENGTH_UNIT=1.496e13
endif

ifeq ($(SPH_SIMULATION),1)
OFLAGS +=-DSPH_TREE
OBJ += getparam.o treegrav.o makecells.o treeio.o clib.o mathfns.o treeload.o

endif


# Code compilation
# ----------------------------------------------------------------------------
%.o: %.c
	$(CC) $(OPT) $(GFLAGS) $(SFLAGS) $(OFLAGS) -c $< 

ifeq ($(DUST_MIX),1)

phaethon :: $(OBJ) phaethon.o
	$(CC) $(OPT) $(GFLAGS) $(SFLAGS) $(OFLAGS) -o phaethon $(OBJ) phaethon.o -lm
	cp phaethon $(EXEC_DIR)
	cp $(OPACITY_DIR)/ptable.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/opacity.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/mopacity.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/ftable.dat $(EXEC_DIR)
	cp $(ISRF_DIR)/isrf.dat $(EXEC_DIR)
	cat Makefile makefiletail.mk>run.info
	mv run.info $(EXEC_DIR)
	cp $(OPACITY2_DIR)/ptable.dat $(EXEC_DIR)/ptable2.dat
	cp $(OPACITY2_DIR)/opacity.dat $(EXEC_DIR)/opacity2.dat
	cp $(OPACITY2_DIR)/mopacity.dat $(EXEC_DIR)/mopacity2.dat
	cp $(OPACITY2_DIR)/ftable.dat $(EXEC_DIR)/ftable2.dat
endif
	
ifeq ($(DUST_MIX),0)
phaethon :: $(OBJ) phaethon.o
	$(CC) $(OPT) $(GFLAGS) $(SFLAGS) $(OFLAGS) -o phaethon $(OBJ) phaethon.o -lm
	cp phaethon $(EXEC_DIR)
	cp $(OPACITY_DIR)/ptable.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/opacity.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/mopacity.dat $(EXEC_DIR)
	cp $(OPACITY_DIR)/ftable.dat $(EXEC_DIR)
	cp $(ISRF_DIR)/isrf.dat $(EXEC_DIR)
	cat Makefile makefiletail.mk>run.info
	mv run.info $(EXEC_DIR)
endif
	
	
mopacity :: $(OBJ_OPA) mopacity.o
	$(CC) $(OPT) -o mopacity $(OBJ_OPA) mopacity.o -lm
	mv mopacity $(OPACITY_DIR)

pemit :: $(OBJ_OPA) pemit.o
	$(CC) $(OPT) -o pemit $(OBJ_OPA) pemit.o -lm
	mv pemit $(OPACITY_DIR)    	

ptest :: $(OBJ_OPA) ptest.o
	$(CC) $(OPT) -o ptest $(OBJ_OPA) ptest.o -lm
	mv ptest $(OPACITY_DIR)

mopatable :: opacity.o make_opacity_table.o inout.o 
	$(CC) $(OPT) -o make_opacity_table  opacity.o inout.o make_opacity_table.o -lm
	mv make_opacity_table  $(OPACITY_DIR)
 
rprofile :: rprofile.o inout.o grid.params.o params.o opacity.o
	
	$(CC) $(OPT) -o rprofile rprofile.o inout.o  grid.params.o params.o opacity.o -lm 
	
	mv rprofile $(EXEC_DIR)
		

clean:
	\rm *.o





