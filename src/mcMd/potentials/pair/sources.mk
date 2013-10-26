mcMd_potentials_pair_=\
    mcMd/potentials/pair/McPairPotential.cpp \
    mcMd/potentials/pair/MdPairPotential.cpp \
    mcMd/potentials/pair/PairFactory.cpp

mcMd_potentials_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_pair_))
mcMd_potentials_pair_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_potentials_pair_:.cpp=.o))

