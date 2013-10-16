inter_angle_=\
    inter/angle/CosineAngle.cpp \
    inter/angle/CosineSqAngle.cpp \
    inter/angle/HarmonicAngle.cpp \

inter_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_angle_))
inter_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_angle_:.cpp=.o))

