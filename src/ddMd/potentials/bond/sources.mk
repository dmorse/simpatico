ddMd_potentials_bond_=\
    ddMd/potentials/bond/BondPotential.cpp \
    ddMd/potentials/bond/BondFactory.cpp 

ddMd_potentials_bond_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_bond_))
ddMd_potentials_bond_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_potentials_bond_:.cpp=.o))

