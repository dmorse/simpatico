mcMd_trajectoryIos_=\
    mcMd/trajectoryIos/TrajectoryIo.cpp \
    mcMd/trajectoryIos/TrajectoryIoFactory.cpp \
    mcMd/trajectoryIos/DCDTrajectoryIo.cpp \
    mcMd/trajectoryIos/LammpsDumpIo.cpp 

mcMd_trajectoryIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_trajectoryIos_))
mcMd_trajectoryIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_trajectoryIos_:.cpp=.o))
