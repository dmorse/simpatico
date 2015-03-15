tools_analyzers_=\
     tools/analyzers/Analyzer.cpp \
     tools/analyzers/AnalyzerManager.cpp \
     tools/analyzers/AtomMSD.cpp \
     tools/analyzers/TrajectoryWriter.cpp \
     tools/analyzers/LammpsDumpWriter.cpp \
     tools/analyzers/PairEnergy.cpp

tools_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_analyzers_))
tools_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_analyzers_:.cpp=.o))

