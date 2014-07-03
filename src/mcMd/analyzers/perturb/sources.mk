mcMd_analyzers_perturb_=\
    mcMd/analyzers/perturb/PerturbDerivative.cpp 

ifdef UTIL_MPI
mcMd_analyzers_perturb_+=\
    mcMd/analyzers/perturb/BennettsMethod.cpp 
endif

mcMd_analyzers_perturb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_perturb_))
mcMd_analyzers_perturb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_perturb_:.cpp=.o))

