mcMd_potentials_coulomb_=\
   mcMd/potentials/coulomb/CoulombPotential.cpp  \
   mcMd/potentials/coulomb/CoulombSystemMixIn.cpp  \
   mcMd/potentials/coulomb/EwaldCoulombPair.cpp \
   mcMd/potentials/coulomb/EwaldCoulombPotential.cpp  \
   mcMd/potentials/coulomb/FineEwaldCoulombPotential.cpp  \
   mcMd/potentials/coulomb/RefinedEwaldCoulombPotential.cpp  \
   mcMd/potentials/coulomb/CoulombFactory.cpp 

mcMd_potentials_coulomb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_coulomb_))
mcMd_potentials_coulomb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_coulomb_:.cpp=.o))

