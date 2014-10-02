spAn_trajectory_=\
   spAn/trajectory/TrajectoryReader.cpp \
   spAn/trajectory/LammpsDumpReader.cpp \
   spAn/trajectory/DdMdTrajectoryReader.cpp \
   spAn/trajectory/TrajectoryReaderFactory.cpp 

spAn_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_trajectory_))
spAn_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_trajectory_:.cpp=.o))

