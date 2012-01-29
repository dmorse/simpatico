
mcMd_mcSimulation_SRCS=$(SRC_DIR)/mcMd/mcSimulation/McDiagnosticManager.cpp \
    $(SRC_DIR)/mcMd/mcSimulation/McModule.cpp \
    $(SRC_DIR)/mcMd/mcSimulation/McSimulation.cpp \
    $(SRC_DIR)/mcMd/mcSimulation/McSystem.cpp 

mcMd_mcSimulation_OBJS=$(mcMd_mcSimulation_SRCS:.cpp=.o)

