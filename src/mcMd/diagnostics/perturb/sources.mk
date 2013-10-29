mcMd_diagnostics_perturb_=\
    mcMd/diagnostics/perturb/PerturbDerivative.cpp 

ifdef UTIL_MPI
mcMd_diagnostics_perturb_+=\
    mcMd/diagnostics/perturb/BennettsMethod.cpp 
endif

mcMd_diagnostics_perturb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_perturb_))
mcMd_diagnostics_perturb_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_diagnostics_perturb_:.cpp=.o))

