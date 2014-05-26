mdCf_tests_configIos_=mdCf/tests/configIos/Test.cc

mdCf_tests_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_tests_configIos_))
mdCf_tests_configIos_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_tests_configIos_:.cc=.o))

