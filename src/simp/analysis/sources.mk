
simp_analysis_= \
    simp/analysis/AnalyzerMixIn.cpp \
    simp/analysis/AverageMixIn.cpp \
    simp/analysis/AverageListMixIn.cpp 

simp_analysis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_analysis_))
simp_analysis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_analysis_:.cpp=.o))

