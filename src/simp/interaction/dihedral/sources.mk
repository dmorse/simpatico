simp_interaction_dihedral_=\
    simp/interaction/dihedral/CosineDihedral.cpp \
    simp/interaction/dihedral/FourierDihedral.cpp 

simp_interaction_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_dihedral_))
simp_interaction_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_dihedral_:.cpp=.o))

