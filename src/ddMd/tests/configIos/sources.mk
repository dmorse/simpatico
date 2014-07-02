ddMd_tests_configIos_=\
   ddMd/tests/configIos/ConfigIoTest.cc \
   ddMd/tests/configIos/SerializeConfigIoTest.cc

ddMd_tests_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_configIos_))
ddMd_tests_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_configIos_:.cc=.o))

