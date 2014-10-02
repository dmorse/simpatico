spAn_analyzers_=\
     spAn/analyzers/Analyzer.cpp \
     spAn/analyzers/AnalyzerFactory.cpp \
     spAn/analyzers/AnalyzerManager.cpp \
     spAn/analyzers/AtomMSD.cpp 

spAn_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_analyzers_))
spAn_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_analyzers_:.cpp=.o))

