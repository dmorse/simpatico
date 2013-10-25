modules_hoomd_potentials_pair_=\
    modules/hoomd/potentials/pair/HoomdPairFactory.cpp \
    modules/hoomd/potentials/pair/HoomdLJPair.cpp \
    modules/hoomd/potentials/pair/HoomdDpdPair.cpp \
    modules/hoomd/potentials/pair/HoomdLJShiftedForcePair.cpp 


modules_hoomd_potentials_pair_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_pair_))
modules_hoomd_potentials_pair_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_pair_:.cpp=.o))

#modules_hoomd_potentials_pair_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_pair_NVCC_))
#modules_hoomd_potentials_pair_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_pair_NVCC_:.cpp=.o))

