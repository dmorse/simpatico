mcMd_diagnostics_simulation_=\
    mcMd/diagnostics/simulation/LogProgress.cpp 

mcMd_diagnostics_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_diagnostics_simulation_))
mcMd_diagnostics_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_diagnostics_simulation_:.cpp=.o))

