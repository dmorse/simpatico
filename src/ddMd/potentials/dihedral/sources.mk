ddMd_potentials_dihedral_=\
    ddMd/potentials/dihedral/DihedralPotential.cpp \
    ddMd/potentials/dihedral/DihedralFactory.cpp 

ddMd_potentials_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_dihedral_))
ddMd_potentials_dihedral_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_potentials_dihedral_:.cpp=.o))

