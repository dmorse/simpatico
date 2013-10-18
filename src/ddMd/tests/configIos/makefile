BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/ddMd/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/configIos/sources.mk

all: $(ddMd_tests_configIos_OBJS) ConfigIoTest SerializeConfigIoTest

clean:
	rm -f $(ddMd_tests_configIos_OBJS) $(ddMd_tests_configIos_OBJS:.o=.d)
	rm -f ConfigIoTest SerializeConfigIoTest

clean-deps:
	rm -f $(ddMd_tests_configIos_OBJS:.o=.d)

ConfigIoTest: ConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o ConfigIoTest ConfigIoTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

SerializeConfigIoTest: SerializeConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o SerializeConfigIoTest SerializeConfigIoTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_configIos_OBJS:.o=.d)

