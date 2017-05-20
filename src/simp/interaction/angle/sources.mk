simp_interaction_angle_=\
    simp/interaction/angle/CosineAngle.cpp \
    simp/interaction/angle/CosineSqAngle.cpp \
    simp/interaction/angle/HarmonicAngle.cpp \

simp_interaction_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_interaction_angle_))
simp_interaction_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_interaction_angle_:.cpp=.o))

