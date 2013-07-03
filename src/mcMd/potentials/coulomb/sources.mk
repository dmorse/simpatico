mcMd_potentials_coulomb_SRCS=\
    $(SRC_DIR)/mcMd/potentials/coulomb/CoulombPotential.cpp  \
    $(SRC_DIR)/mcMd/potentials/coulomb/EwaldCoulombPair.cpp  \
    $(SRC_DIR)/mcMd/potentials/coulomb/MdCoulombPotential.cpp  

mcMd_potentials_coulomb_OBJS=$(mcMd_potentials_coulomb_SRCS:.cpp=.o)

