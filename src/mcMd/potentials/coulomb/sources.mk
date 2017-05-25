mcMd_potentials_coulomb_=\
   mcMd/potentials/coulomb/CoulombFactory.cpp  \
   mcMd/potentials/coulomb/MdCoulombPotential.cpp \
   mcMd/potentials/coulomb/MdEwaldPotential.cpp \
   mcMd/potentials/coulomb/EwaldRSpaceAccumulator.cpp 

ifdef SIMP_FFTW
mcMd_potentials_coulomb_+=\
   mcMd/potentials/coulomb/MdSpmePotential.cpp 
endif

mcMd_potentials_coulomb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_coulomb_))
mcMd_potentials_coulomb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_coulomb_:.cpp=.o))

