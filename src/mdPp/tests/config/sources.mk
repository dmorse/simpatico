mdPp_tests_config_=mdPp/tests/config/Test.cc

mdPp_tests_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_config_))
mdPp_tests_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_config_:.cc=.o))

