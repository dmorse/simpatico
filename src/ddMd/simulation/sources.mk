ddMd_simulation_= \
     ddMd/simulation/Simulation.cpp \
     ddMd/simulation/SimulationAccess.cpp 

ddMd_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_simulation_))
ddMd_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_simulation_:.cpp=.o))

