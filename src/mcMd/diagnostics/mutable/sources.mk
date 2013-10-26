mcMd_diagnostics_mutable_=\
    mcMd/diagnostics/mutable/SemiGrandDistribution.cpp \
    mcMd/diagnostics/mutable/TypeDistribution.cpp 

mcMd_diagnostics_mutable_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_mutable_))
mcMd_diagnostics_mutable_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_diagnostics_mutable_:.cpp=.o))

