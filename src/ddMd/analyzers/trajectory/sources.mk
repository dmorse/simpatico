ddMd_analyzers_trajectory_=\
     ddMd/analyzers/trajectory/ConfigWriter.cpp\
     ddMd/analyzers/trajectory/TrajectoryWriter.cpp\
     ddMd/analyzers/trajectory/DdMdTrajectoryWriter.cpp\
     ddMd/analyzers/trajectory/DdMdGroupTrajectoryWriter.cpp\
     ddMd/analyzers/trajectory/LammpsDumpWriter.cpp

ddMd_analyzers_trajectory_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_trajectory_))
ddMd_analyzers_trajectory_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_trajectory_:.cpp=.o))

