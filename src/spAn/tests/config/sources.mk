spAn_tests_config_=spAn/tests/config/Test.cc

spAn_tests_config_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_config_))
spAn_tests_config_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_config_:.cc=.o))

