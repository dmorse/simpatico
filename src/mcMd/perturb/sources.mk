include $(SRC_DIR)/mcMd/perturb/mcSystem/sources.mk

mcMd_perturb_SRCS=$(mcMd_perturb_mcSystem_SRCS) \
    $(SRC_DIR)/mcMd/perturb/Perturbation.cpp \
    $(SRC_DIR)/mcMd/perturb/ReplicaMove.cpp 

mcMd_perturb_OBJS=$(mcMd_perturb_SRCS:.cpp=.o)

