mdPp_analyzers_=\
     mdPp/analyzers/Analyzer.cpp \
     mdPp/analyzers/AnalyzerFactory.cpp \
     mdPp/analyzers/AnalyzerManager.cpp

mdPp_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_analyzers_))
mdPp_analyzers_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_analyzers_:.cpp=.o))

