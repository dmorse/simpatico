modules_hoomd_potentials_bond_=

modules_hoomd_potentials_bond_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_bond_))
modules_hoomd_potentials_bond_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_bond_:.cpp=.o))

#modules_hoomd_potentials_bond_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_potentials_bond_NVCC_))
#modules_hoomd_potentials_bond_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_potentials_bond_NVCC_:.cpp=.o))

