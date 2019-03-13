ddMd_analyzers_base_=\
     ddMd/analyzers/base/AverageAnalyzer.cpp\
     ddMd/analyzers/base/TensorAverageAnalyzer.cpp\
     ddMd/analyzers/base/SymmTensorAverageAnalyzer.cpp\
     ddMd/analyzers/base/AverageListAnalyzer.cpp

ddMd_analyzers_base_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_base_))
ddMd_analyzers_base_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_base_:.cpp=.o))

