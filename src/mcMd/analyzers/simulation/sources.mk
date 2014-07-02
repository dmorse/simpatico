mcMd_analyzers_simulation_=\
    mcMd/analyzers/simulation/LogProgress.cpp 

mcMd_analyzers_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_simulation_))
mcMd_analyzers_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_simulation_:.cpp=.o))

