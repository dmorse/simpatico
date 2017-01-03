include $(SRC_DIR)/modules/hoomd/analyzers/sources.mk
include $(SRC_DIR)/modules/hoomd/mcMoves/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/pair/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/bond/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/link/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/external/sources.mk
include $(SRC_DIR)/modules/hoomd/perturbation/sources.mk

#-----------------------------------------------------------------
# CPP files

modules_hoomd_= \
    $(modules_hoomd_analyzers_) \
    $(modules_hoomd_mcMoves_) \
    $(modules_hoomd_potentials_pair_) \
    $(modules_hoomd_potentials_bond_) \
    $(modules_hoomd_potentials_link_) \
    $(modules_hoomd_potentials_external_) \
    $(modules_hoomd_perturbation_) \
    modules/hoomd/HoomdMcModule.cpp 

modules_hoomd_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_))
modules_hoomd_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_:.cpp=.o))

#-----------------------------------------------------------------
# Cuda C files

modules_hoomd_NVCC_=\
    $(modules_hoomd_analyzers_NVCC_)

modules_hoomd_NVCC_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_NVCC_))
modules_hoomd_NVCC_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_NVCC_:.cu=.cu.o))

