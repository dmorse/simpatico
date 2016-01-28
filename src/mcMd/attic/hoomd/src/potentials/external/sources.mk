modules_hoomd_potentials_external_=\
    modules/hoomd/potentials/external/HoomdExternalFactory.cpp \
    modules/hoomd/potentials/external/HoomdPeriodicExternal.cpp \
    modules/hoomd/potentials/external/HoomdLocalLamellarOrderingExternal.cpp 

modules_hoomd_potentials_external_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_external_))
modules_hoomd_potentials_external_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_external_:.cpp=.o))

#modules_hoomd_potentials_external_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_external_NVCC_))
#modules_hoomd_potentials_external_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_external_NVCC_:.cpp=.o))

