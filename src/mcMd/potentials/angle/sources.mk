mcMd_potentials_angle_=\
    mcMd/potentials/angle/AnglePotential.cpp \
    mcMd/potentials/angle/AngleFactory.cpp

mcMd_potentials_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_potentials_angle_))
mcMd_potentials_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_potentials_angle_:.cpp=.o))

