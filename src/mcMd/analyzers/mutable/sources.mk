mcMd_analyzers_mutable_=\
    mcMd/analyzers/mutable/SemiGrandDistribution.cpp \
    mcMd/analyzers/mutable/TypeDistribution.cpp 

mcMd_analyzers_mutable_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_mutable_))
mcMd_analyzers_mutable_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_mutable_:.cpp=.o))

