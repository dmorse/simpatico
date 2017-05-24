simp_interaction_coulomb_=\
    simp/interaction/coulomb/EwaldInteraction.cpp

simp_interaction_coulomb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_coulomb_))
simp_interaction_coulomb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_coulomb_:.cpp=.o))

