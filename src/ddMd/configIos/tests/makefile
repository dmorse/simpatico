SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include sources.mk

all: $(ddMd_configIos_tests_OBJS) ConfigIoTest

clean:
	rm -f $(ddMd_configIos_tests_OBJS) $(ddMd_configIos_tests_OBJS:.o=.d)
	rm -f ConfigIoTest

clean-deps:
	rm -f $(ddMd_configIos_tests_OBJS:.o=.d)

ConfigIoTest: ConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(UTIL_DEFS) $(INTER_DEFS) $(DDMD_DEFS) \
               -o ConfigIoTest ConfigIoTest.o $(ddMd_LIB) $(util_LIB)

-include $(ddMd_configIos_tests_OBJS:.o=.d)

