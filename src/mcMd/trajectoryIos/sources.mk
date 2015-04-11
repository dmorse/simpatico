mcMd_trajectory_=\
    mcMd/trajectoryIos/TrajectoryReader.cpp \
    mcMd/trajectoryIos/TrajectoryReaderFactory.cpp \
    mcMd/trajectoryIos/DCDTrajectoryReader.cpp \
    mcMd/trajectoryIos/LammpsDumpReader.cpp \
    mcMd/trajectoryIos/DdMdTrajectoryReader.cpp 

mcMd_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_trajectory_))
mcMd_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_trajectory_:.cpp=.o))

