mcMd_analyzers_trajectory_=\
    mcMd/analyzers/trajectory/TrajectoryWriter.cpp \
    mcMd/analyzers/trajectory/ConfigDumpWriter.cpp \
    mcMd/analyzers/trajectory/DdMdTrajectoryWriter.cpp \
    mcMd/analyzers/trajectory/SmpTrajectoryWriter.cpp \

mcMd_analyzers_trajectory_SRCS=\
    $(addprefix $(SRC_DIR)/, $(mcMd_analyzers_trajectory_))
mcMd_analyzers_trajectory_OBJS=\
    $(addprefix $(BLD_DIR)/, $(mcMd_analyzers_trajectory_:.cpp=.o))

