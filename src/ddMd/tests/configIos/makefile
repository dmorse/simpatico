SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/ddMd/defines.mk
include $(TESTS_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk


include sources.mk

all: $(tests_ddMd_configIos_OBJS) ConfigIoTest SerializeConfigIoTest

clean:
	rm -f $(tests_ddMd_configIos_OBJS) $(tests_ddMd_configIos_OBJS:.o=.d)
	rm -f ConfigIoTest SerializeConfigIoTest
	rm -f out

clean-deps:
	rm -f $(tests_ddMd_configIos_OBJS:.o=.d)

ConfigIoTest: ConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ConfigIoTest ConfigIoTest.o $(ddMd_LIB) $(util_LIB)

SerializeConfigIoTest: SerializeConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o SerializeConfigIoTest SerializeConfigIoTest.o $(ddMd_LIB) $(util_LIB)

-include $(tests_ddMd_configIos_OBJS:.o=.d)

