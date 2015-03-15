tools_trajectory_=\
   tools/trajectory/TrajectoryReader.cpp \
   tools/trajectory/LammpsDumpReader.cpp \
   tools/trajectory/DdMdTrajectoryReader.cpp \
   tools/trajectory/TrajectoryReaderFactory.cpp 

tools_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_trajectory_))
tools_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_trajectory_:.cpp=.o))

