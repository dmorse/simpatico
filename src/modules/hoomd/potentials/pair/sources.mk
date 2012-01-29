modules_hoomd_potentials_pair_SRCS=$(SRC_DIR)/modules/hoomd/potentials/pair/HoomdPairFactory.cpp \
    $(SRC_DIR)/modules/hoomd/potentials/pair/HoomdLJPair.cpp \
    $(SRC_DIR)/modules/hoomd/potentials/pair/HoomdDpdPair.cpp \
    $(SRC_DIR)/modules/hoomd/potentials/pair/HoomdLJShiftedForcePair.cpp \


modules_hoomd_potentials_pair_OBJS=$(modules_hoomd_potentials_pair_SRCS:.cpp=.o)

