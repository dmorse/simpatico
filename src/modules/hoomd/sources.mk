include $(SRC_DIR)/modules/hoomd/diagnostics/sources.mk
include $(SRC_DIR)/modules/hoomd/mcMoves/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/pair/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/bond/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/link/sources.mk
include $(SRC_DIR)/modules/hoomd/potentials/external/sources.mk
include $(SRC_DIR)/modules/hoomd/perturbation/sources.mk

modules_hoomd_SRCS= \
    $(modules_hoomd_diagnostics_SRCS) \
    $(modules_hoomd_mcMoves_SRCS) \
    $(modules_hoomd_potentials_pair_SRCS) \
    $(modules_hoomd_potentials_bond_SRCS) \
    $(modules_hoomd_potentials_link_SRCS) \
    $(modules_hoomd_potentials_external_SRCS) \
    $(modules_hoomd_perturbation_SRCS) \
    $(SRC_DIR)/modules/hoomd/HoomdMcModule.cpp 

modules_hoomd_OBJS=$(modules_hoomd_SRCS:.cpp=.o)

modules_hoomd_NVCC_SRCS = $(modules_hoomd_diagnostics_NVCC_SRCS)

modules_hoomd_NVCC_OBJS = $(modules_hoomd_NVCC_SRCS:.cu=.cu.o)
