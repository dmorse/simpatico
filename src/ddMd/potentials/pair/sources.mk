ddMd_potentials_pair_= \
    ddMd/potentials/pair/PairPotential.cpp \
    ddMd/potentials/pair/PairFactory.cpp 

ddMd_potentials_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_pair_))
ddMd_potentials_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_potentials_pair_:.cpp=.o))

