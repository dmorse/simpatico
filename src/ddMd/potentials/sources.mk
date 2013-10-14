include $(SRC_DIR)/ddMd/potentials/pair/sources.mk

ddMd_potentials_SRCS=\
    $(SRC_DIR)/ddMd/potentials/Potential.cpp\
    $(ddMd_potentials_pair_SRCS) 

ifdef INTER_BOND
include $(SRC_DIR)/ddMd/potentials/bond/sources.mk
ddMd_potentials_SRCS+=$(ddMd_potentials_bond_SRCS)
endif

ifdef INTER_ANGLE
include $(SRC_DIR)/ddMd/potentials/angle/sources.mk
ddMd_potentials_SRCS+=$(ddMd_potentials_angle_SRCS)
endif

ifdef INTER_DIHEDRAL
include $(SRC_DIR)/ddMd/potentials/dihedral/sources.mk
ddMd_potentials_SRCS+=$(ddMd_potentials_dihedral_SRCS)
endif

ifdef INTER_EXTERNAL
include $(SRC_DIR)/ddMd/potentials/external/sources.mk
ddMd_potentials_SRCS+=$(ddMd_potentials_external_SRCS)
endif

ddMd_potentials_OBJS=$(ddMd_potentials_SRCS:.cpp=.o)
