
mcMd_simulation_SRCS=$(SRC_DIR)/mcMd/simulation/McMd_mpi.cpp \
    $(SRC_DIR)/mcMd/simulation/Simulation.cpp \
    $(SRC_DIR)/mcMd/simulation/SubSystem.cpp \
    $(SRC_DIR)/mcMd/simulation/System.cpp 

mcMd_simulation_OBJS=$(mcMd_simulation_SRCS:.cpp=.o)

