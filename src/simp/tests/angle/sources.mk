simp_tests_angle_=simp/tests/angle/Test.cc


simp_tests_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_angle_))
simp_tests_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_angle_:.cc=.o))

