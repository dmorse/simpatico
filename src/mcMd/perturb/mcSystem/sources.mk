mcMd_perturb_mcSystem_SRCS=\
    $(SRC_DIR)/mcMd/perturb/mcSystem/McEnergyPerturbation.cpp \
    $(SRC_DIR)/mcMd/perturb/mcSystem/McPerturbationFactory.cpp 

ifdef INTER_EXTERNAL
mcMd_perturb_mcSystem_SRCS+=\
    $(SRC_DIR)/mcMd/perturb/mcSystem/McExternalPerturbation.cpp 
endif

mcMd_perturb_mcSystem_OBJS=$(mcMd_perturb_mcSystem_SRCS:.cpp=.o)

