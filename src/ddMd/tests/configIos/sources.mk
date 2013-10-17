tests_ddMd_configIos_SRCS=\
   $(TESTS_DIR)/ddMd/configIos/ConfigIoTest.cc \
   $(TESTS_DIR)/ddMd/configIos/SerializeConfigIoTest.cc

tests_ddMd_configIos_OBJS=$(tests_ddMd_configIos_SRCS:.cc=.o)

