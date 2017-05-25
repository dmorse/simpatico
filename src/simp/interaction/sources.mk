
simp_interaction_=

ifndef SIMP_NOPAIR
include $(SRC_DIR)/simp/interaction/pair/sources.mk
simp_interaction_+=$(simp_interaction_pair_)
endif

ifdef SIMP_BOND
include $(SRC_DIR)/simp/interaction/bond/sources.mk
simp_interaction_+=$(simp_interaction_bond_)
endif

ifdef SIMP_ANGLE
include $(SRC_DIR)/simp/interaction/angle/sources.mk
simp_interaction_+=$(simp_interaction_angle_)
endif

ifdef SIMP_DIHEDRAL
include $(SRC_DIR)/simp/interaction/dihedral/sources.mk
simp_interaction_+=$(simp_interaction_dihedral_)
endif

ifdef SIMP_COULOMB
include $(SRC_DIR)/simp/interaction/coulomb/sources.mk
simp_interaction_+=$(simp_interaction_coulomb_)
endif

ifdef SIMP_EXTERNAL
include $(SRC_DIR)/simp/interaction/external/sources.mk
simp_interaction_+=$(simp_interaction_external_)
endif

simp_interaction_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_))
simp_interaction_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_:.cpp=.o))


