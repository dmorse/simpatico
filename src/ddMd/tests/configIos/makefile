OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/configIos/sources.mk

all: $(ddMd_tests_configIos_OBJS) \
     ConfigIoTest \
     SerializeConfigIoTest 

clean:
	rm -f $(ddMd_tests_configIos_OBJS) $(ddMd_tests_configIos_OBJS:.o=.d)
	rm -f ConfigIoTest SerializeConfigIoTest
	rm -f binary out

ConfigIoTest: ConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) -o ConfigIoTest ConfigIoTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

SerializeConfigIoTest: SerializeConfigIoTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) -o SerializeConfigIoTest SerializeConfigIoTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_configIos_OBJS:.o=.d)

