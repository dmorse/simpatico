spAn_tests_configIos_=spAn/tests/configIos/Test.cc

spAn_tests_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_configIos_))
spAn_tests_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_configIos_:.cc=.o))

