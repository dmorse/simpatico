mdPp_analyzers_=\
     mdPp/analyzers/Analyzer.cpp \
     mdPp/analyzers/AnalyzerManager.cpp \
     mdPp/analyzers/AtomMSD.cpp \
     mdPp/analyzers/TrajectoryWriter.cpp \
     mdPp/analyzers/LammpsDumpWriter.cpp \
     mdPp/analyzers/PairEnergy.cpp

mdPp_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_analyzers_))
mdPp_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_analyzers_:.cpp=.o))

