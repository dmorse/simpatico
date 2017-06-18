
mcMd_mcSimulation_=\
    mcMd/mcSimulation/McSimulation.cpp \
    mcMd/mcSimulation/McSystem.cpp \
    mcMd/mcSimulation/McSystemInterface.cpp \
    mcMd/mcSimulation/McAnalyzerManager.cpp \
    mcMd/mcSimulation/McCommandManager.cpp \

mcMd_mcSimulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mcSimulation_))
mcMd_mcSimulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mcSimulation_:.cpp=.o))

