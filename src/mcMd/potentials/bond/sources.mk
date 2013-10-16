mcMd_potentials_bond_=\
    mcMd/potentials/bond/BondFactory.cpp \
    mcMd/potentials/bond/BondPotential.cpp 

mcMd_potentials_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_bond_))
mcMd_potentials_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_bond_:.cpp=.o))

