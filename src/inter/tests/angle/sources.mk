inter_tests_angle_=inter/tests/angle/Test.cc


inter_tests_angle_SRCS=\
     $(addprefix $(SRC_DIR)/, $(inter_tests_angle_))
inter_tests_angle_OBJS=\
     $(addprefix $(BLD_DIR)/, $(inter_tests_angle_:.cc=.o))

