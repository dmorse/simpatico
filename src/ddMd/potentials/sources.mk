include $(SRC_DIR)/ddMd/potentials/pair/sources.mk

ddMd_potentials_= \
   ddMd/potentials/Potential.cpp \
   $(ddMd_potentials_pair_) 

ifdef SIMP_BOND
include $(SRC_DIR)/ddMd/potentials/bond/sources.mk
ddMd_potentials_+=$(ddMd_potentials_bond_)
endif

ifdef SIMP_ANGLE
include $(SRC_DIR)/ddMd/potentials/angle/sources.mk
ddMd_potentials_+=$(ddMd_potentials_angle_)
endif

ifdef SIMP_DIHEDRAL
include $(SRC_DIR)/ddMd/potentials/dihedral/sources.mk
ddMd_potentials_+=$(ddMd_potentials_dihedral_)
endif

ifdef SIMP_EXTERNAL
include $(SRC_DIR)/ddMd/potentials/external/sources.mk
ddMd_potentials_+=$(ddMd_potentials_external_)
endif

ddMd_potentials_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_))
ddMd_potentials_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_potentials_:.cpp=.o))

