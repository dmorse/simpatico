mdPp_trajectory_=\
   mdPp/trajectory/TrajectoryReader.cpp \
   mdPp/trajectory/LammpsDumpReader.cpp \
   mdPp/trajectory/DdMdTrajectoryReader.cpp \
   mdPp/trajectory/TrajectoryReaderFactory.cpp 

mdPp_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_trajectory_))
mdPp_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_trajectory_:.cpp=.o))

