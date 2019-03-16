
simp_analyzers_= \
    simp/analyzers/AnalyzerMixIn.cpp \
    simp/analyzers/AverageMixIn.cpp 

simp_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_analyzers_))
simp_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_analyzers_:.cpp=.o))

