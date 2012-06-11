SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(tests_ddMd_configIos_OBJS) ConfigIoTest

clean:
	rm -f $(tests_ddMd_configIos_OBJS) $(tests_ddMd_configIos_OBJS:.o=.d)
	rm -f ConfigIoTest
	rm -f out

clean-deps:
	rm -f $(tests_ddMd_configIos_OBJS:.o=.d)

ConfigIoTest: ConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ConfigIoTest ConfigIoTest.o $(ddMd_LIB) $(util_LIB)

-include $(tests_ddMd_configIos_OBJS:.o=.d)

