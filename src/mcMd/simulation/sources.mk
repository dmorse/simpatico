
mcMd_simulation_= \
    mcMd/simulation/Simulation.cpp \
    mcMd/simulation/System.cpp \
    mcMd/simulation/SystemInterface.cpp \

ifdef UTIL_MPI
mcMd_simulation_ += mcMd/simulation/McMd_mpi.cpp 
endif

mcMd_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_simulation_))
mcMd_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_simulation_:.cpp=.o))

