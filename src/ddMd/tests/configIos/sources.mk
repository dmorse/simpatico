tests_ddMd_configIos_SRCS=\
   $(TESTS_DIR)/ddMd/configIos/ConfigIoTest.cc \
   $(TESTS_DIR)/ddMd/configIos/SerializeConfigIoTest.cc

tests_ddMd_configIos_OBJS=$(tests_ddMd_configIos_SRCS:.cc=.o)

ddMd_tests_TMP_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_TMP_))
ddMd_tests_TMP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_TMP_:.cpp=.o))

