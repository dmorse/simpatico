mcMd_potentials_bond_SRCS=$(SRC_DIR)/mcMd/potentials/bond/BondFactory.cpp \
    $(SRC_DIR)/mcMd/potentials/bond/BondPotential.cpp \
    $(SRC_DIR)/mcMd/potentials/bond/FeneBond.cpp \
    $(SRC_DIR)/mcMd/potentials/bond/HarmonicBond.cpp \
    $(SRC_DIR)/mcMd/potentials/bond/HarmonicL0Bond.cpp 

mcMd_potentials_bond_OBJS=$(mcMd_potentials_bond_SRCS:.cpp=.o)

