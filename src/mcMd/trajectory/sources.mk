mcMd_trajectory_=\
    mcMd/trajectory/TrajectoryReader.cpp \
    mcMd/trajectory/TrajectoryReaderFactory.cpp \
    mcMd/trajectory/DCDTrajectoryReader.cpp \
    mcMd/trajectory/LammpsDumpReader.cpp \
    mcMd/trajectory/DdMdTrajectoryReader.cpp 

mcMd_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_trajectory_))
mcMd_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_trajectory_:.cpp=.o))

