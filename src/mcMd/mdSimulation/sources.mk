
mcMd_mdSimulation_=\
    mcMd/mdSimulation/MdSimulation.cpp \
    mcMd/mdSimulation/MdSystem.cpp \
    mcMd/mdSimulation/MdSystemInterface.cpp \
    mcMd/mdSimulation/MdAnalyzerManager.cpp \
    mcMd/mdSimulation/MdCommandManager.cpp \

mcMd_mdSimulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_mdSimulation_))
mcMd_mdSimulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_mdSimulation_:.cpp=.o))

