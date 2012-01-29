
mcMd_mdSimulation_SRCS=$(SRC_DIR)/mcMd/mdSimulation/MdDiagnosticManager.cpp \
    $(SRC_DIR)/mcMd/mdSimulation/MdModule.cpp \
    $(SRC_DIR)/mcMd/mdSimulation/MdSimulation.cpp \
    $(SRC_DIR)/mcMd/mdSimulation/MdSystem.cpp 

mcMd_mdSimulation_OBJS=$(mcMd_mdSimulation_SRCS:.cpp=.o)

