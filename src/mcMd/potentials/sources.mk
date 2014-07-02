include $(SRC_DIR)/mcMd/potentials/bond/sources.mk
mcMd_potentials_=\
    $(mcMd_potentials_bond_) 

ifndef INTER_NOPAIR
include $(SRC_DIR)/mcMd/potentials/pair/sources.mk
mcMd_potentials_+=$(mcMd_potentials_pair_) 
endif

ifdef INTER_ANGLE
include $(SRC_DIR)/mcMd/potentials/angle/sources.mk
mcMd_potentials_+=$(mcMd_potentials_angle_) 
endif

ifdef INTER_DIHEDRAL
include $(SRC_DIR)/mcMd/potentials/dihedral/sources.mk
mcMd_potentials_+=$(mcMd_potentials_dihedral_) 
endif

ifdef MCMD_LINK
include $(SRC_DIR)/mcMd/potentials/link/sources.mk
mcMd_potentials_+=$(mcMd_potentials_link_) 
endif

ifdef INTER_EXTERNAL
include $(SRC_DIR)/mcMd/potentials/external/sources.mk
mcMd_potentials_+=$(mcMd_potentials_external_) 
endif

mcMd_potentials_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_))
mcMd_potentials_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_:.cpp=.o))

