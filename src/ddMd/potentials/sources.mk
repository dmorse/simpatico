include $(SRC_DIR)/ddMd/potentials/bond/sources.mk
include $(SRC_DIR)/ddMd/potentials/pair/sources.mk

ddMd_potentials_SRCS=$(ddMd_potentials_bond_SRCS) \
    $(ddMd_potentials_pair_SRCS) \
    $(SRC_DIR)/ddMd/potentials/PairPotential.cpp \
    $(SRC_DIR)/ddMd/potentials/BondPotential.cpp 

ddMd_potentials_OBJS=$(ddMd_potentials_SRCS:.cpp=.o)

