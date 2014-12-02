tools_tests_config_=tools/tests/config/Test.cc

tools_tests_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_config_))
tools_tests_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_config_:.cc=.o))

