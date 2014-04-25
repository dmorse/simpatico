mdPp_tests_configIos_=mdPp/tests/configIos/Test.cc

mdPp_tests_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_configIos_))
mdPp_tests_configIos_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_tests_configIos_:.cc=.o))

