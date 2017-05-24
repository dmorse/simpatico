simp_tests_interaction_angle_=\
     simp/tests/interaction/angle/Test.cc

simp_tests_interaction_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_interaction_angle_))
simp_tests_interaction_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_interaction_angle_:.cc=.o))

