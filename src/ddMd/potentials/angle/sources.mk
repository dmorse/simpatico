ddMd_potentials_angle_=\
    ddMd/potentials/angle/AnglePotential.cpp \
    ddMd/potentials/angle/AngleFactory.cpp 

ddMd_potentials_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_potentials_angle_))
ddMd_potentials_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_potentials_angle_:.cpp=.o))

