mcMd_analyzers_base_= \
     mcMd/analyzers/base/SystemAnalyzer.cpp \
     mcMd/analyzers/base/AverageAnalyzer.cpp \
     mcMd/analyzers/base/AverageListAnalyzer.cpp \
     mcMd/analyzers/base/DistributionAnalyzer.cpp 

mcMd_analyzers_base_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_base_))
mcMd_analyzers_base_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_base_:.cpp=.o))

