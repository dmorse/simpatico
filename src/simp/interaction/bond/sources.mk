simp_interaction_bond_=\
    simp/interaction/bond/FeneBond.cpp \
    simp/interaction/bond/HarmonicBond.cpp \
    simp/interaction/bond/HarmonicL0Bond.cpp 

simp_interaction_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_bond_))
simp_interaction_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_bond_:.cpp=.o))

