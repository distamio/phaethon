# -----------------------------------------------------------------------------
# PHAETHON MAKEFILE 
# -----------------------------------------------------------------------------

CC                		=gcc 
OPENMP             		=0
FASTMATH			=1
OPTIMISE			=3
PROFILE                        	=0 
MAIN_DIR             		=$(PWD)/main
GRID_DIR           		=$(PWD)/grid.1
#OPACITY_DIR                    =$(PWD)/opacity.disc.wood02_0.2_1200
#OPACITY_DIR                    =$(PWD)/opacity.OH5_0.2_1200
 OPACITY_DIR                    =$(PWD)/opacity.OH5_0.1_100
#OPACITY_DIR                    =$(PWD)/opacity.OH8_0.2_1200
#OPACITY_DIR                    =$(PWD)/opacity.DRAINE03.MRN.0.2_1200
#OPACITY2_DIR                   =$(PWD)/opacity.disc.wood02_0.2_1200
OPACITY_UTILS_DIR               =$(PWD)/opacity.utils
UTILS_DIR                       =$(PWD)/utils
ISRF_DIR                        =$(PWD)/isrf_andre
EXEC_DIR                        =$(PWD)/execdir

# -----------------------------------------------------------------------------
# GRID 
# -----------------------------------------------------------------------------
SPHERICAL          	 	=1
DIM		           	=2
# -----------------------------------------------------------------------------
# SEDs  & INTENSITY MAPS
# ----------------------------------------------------------------------------
NFBINS				=200  # 200 #400  #330 120
FS				=5
NLBINS				=1
NMBINS				=1 

DO_ISOMAP	   		=1
NBBINS				=80
NXPIXELS 			=128
NYPIXELS 			=128	

ISODIM		    		=2
EXTRA_SYMMETRY 			=0
SPHERICAL_ISOMAP  		=1
SPHERICAL_SED			=1

# -----------------------------------------------------------------------------
# DUST PROPERTIES 
# -----------------------------------------------------------------------------

DUST_MIX                        = 0   
DUST_KIND			= 0

# -----------------------------------------------------------------------------
# OTHER OPTIONS 
# -----------------------------------------------------------------------------
NO_INTERACTION     		=0 
MULTIPLE_PHOTONS   		=0
TREE_INFO_OUT      		=0
INFO_OUTPUT        		=0 
INC_SHADOW_PHOTONS		=0
# LENGTH_UNIT (PC or AU) for input params.dat file and params.grid.dat
LENGTH_UNIT			=PC 
SPH_SIMULATION          	=0


#################################################################################

# DUST_KIND                     1: DISC_IN_CLOUD


# Now include the Makefile tail which contains all the processed options 
# and generates list of object files for compilation
# ----------------------------------------------------------------------------
include makefiletail.mk

