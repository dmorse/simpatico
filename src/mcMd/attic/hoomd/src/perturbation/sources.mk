modules_hoomd_perturbation_=\
   modules/hoomd/perturbation/HoomdMcPerturbationFactory.cpp 

modules_hoomd_perturbation_SRCS=\
    $(addprefix $(SRC_DIR)/, $(modules_hoomd_perturbation_))
modules_hoomd_perturbation_OBJS=\
    $(addprefix $(BLD_DIR)/, $(modules_hoomd_perturbation_:.cpp=.o))

#modules_hoomd_perturbation_NVCC_SRCS=\
#    $(addprefix $(SRC_DIR)/, $(modules_hoomd_perturbation_NVCC_))
#modules_hoomd_perturbation_NVCC_OBJS=\
#    $(addprefix $(BLD_DIR)/, $(modules_hoomd_perturbation_NVCC_:.cpp=.o))

