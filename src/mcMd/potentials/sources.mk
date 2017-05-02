mcMd_potentials_= 

include $(SRC_DIR)/mcMd/potentials/misc/sources.mk
mcMd_potentials_+=$(mcMd_potentials_misc_) 

ifndef SIMP_NOPAIR
include $(SRC_DIR)/mcMd/potentials/pair/sources.mk
mcMd_potentials_+=$(mcMd_potentials_pair_) 
endif

ifdef SIMP_BOND
include $(SRC_DIR)/mcMd/potentials/bond/sources.mk
mcMd_potentials_+=\
    $(mcMd_potentials_bond_) 
endif

ifdef SIMP_ANGLE
include $(SRC_DIR)/mcMd/potentials/angle/sources.mk
mcMd_potentials_+=$(mcMd_potentials_angle_) 
endif

ifdef SIMP_DIHEDRAL
include $(SRC_DIR)/mcMd/potentials/dihedral/sources.mk
mcMd_potentials_+=$(mcMd_potentials_dihedral_) 
endif

ifdef SIMP_COULOMB
include $(SRC_DIR)/mcMd/potentials/coulomb/sources.mk
mcMd_potentials_+=$(mcMd_potentials_coulomb_) 
endif

ifdef MCMD_LINK
include $(SRC_DIR)/mcMd/potentials/link/sources.mk
mcMd_potentials_+=$(mcMd_potentials_link_) 
endif

ifdef SIMP_EXTERNAL
include $(SRC_DIR)/mcMd/potentials/external/sources.mk
mcMd_potentials_+=$(mcMd_potentials_external_) 
endif

mcMd_potentials_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_))
mcMd_potentials_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_:.cpp=.o))

