mcMd_potentials_dihedral_=\
    mcMd/potentials/dihedral/DihedralPotential.cpp \
    mcMd/potentials/dihedral/DihedralFactory.cpp 

mcMd_potentials_dihedral_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_dihedral_))
mcMd_potentials_dihedral_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_dihedral_:.cpp=.o))

