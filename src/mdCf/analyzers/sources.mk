mdCf_analyzers_=\
     mdCf/analyzers/Analyzer.cpp \
     mdCf/analyzers/AnalyzerFactory.cpp \
     mdCf/analyzers/AnalyzerManager.cpp \
     mdCf/analyzers/AtomMSD.cpp 

mdCf_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_analyzers_))
mdCf_analyzers_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_analyzers_:.cpp=.o))

