ddMd_potentials_bond_SRCS=\
    $(SRC_DIR)/ddMd/potentials/bond/BondPotential.cpp \
    $(SRC_DIR)/ddMd/potentials/bond/BondFactory.cpp \
    $(SRC_DIR)/ddMd/potentials/bond/HarmonicBond.cpp \
    $(SRC_DIR)/ddMd/potentials/bond/HarmonicL0Bond.cpp 

ddMd_potentials_bond_OBJS=$(ddMd_potentials_bond_SRCS:.cpp=.o)

