include $(SRC_DIR)/mcMd/perturb/mcSystem/sources.mk

ifdef MCMD_PERTURB
mcMd_perturb_=\
    $(mcMd_perturb_mcSystem_) \
    mcMd/perturb/Perturbation.cpp 
ifdef UTIL_MPI
mcMd_perturb_+=\
    mcMd/perturb/ReplicaMove.cpp 
endif
endif

mcMd_perturb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_perturb_))
mcMd_perturb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_perturb_:.cpp=.o))

