mcMd_diagnostics_perturb_SRCS=\
    $(SRC_DIR)/mcMd/diagnostics/perturb/PerturbDerivative.cpp 

ifdef UTIL_MPI
mcMd_diagnostics_perturb_SRCS+=\
    $(SRC_DIR)/mcMd/diagnostics/perturb/BennettsMethod.cpp 
endif

mcMd_diagnostics_perturb_OBJS=$(mcMd_diagnostics_perturb_SRCS:.cpp=.o)

