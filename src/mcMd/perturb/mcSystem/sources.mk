mcMd_perturb_mcSystem_=\
    mcMd/perturb/mcSystem/McEnergyPerturbation.cpp \
    mcMd/perturb/mcSystem/McPerturbationFactory.cpp 

mcMd_perturb_mcSystem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_perturb_mcSystem_))
mcMd_perturb_mcSystem_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_perturb_mcSystem_:.cpp=.o))

