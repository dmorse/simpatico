mcMd_analyzers_base_=\
     mcMd/analyzers/base/McAverageAnalyzer.cpp \
     mcMd/analyzers/base/McAverageListAnalyzer.cpp

mcMd_analyzers_base_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_base_))
mcMd_analyzers_base_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_base_:.cpp=.o))

