modules_hoomd_potentials_link_=

modules_hoomd_potentials_link_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_link_))
modules_hoomd_potentials_link_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_link_:.cpp=.o))

#modules_hoomd_potentials_link_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_link_NVCC_))
#modules_hoomd_potentials_link_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_link_NVCC_:.cpp=.o))

