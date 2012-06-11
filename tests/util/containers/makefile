SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
all: $(tests_util_containers_OBJS) Test

clean:
	rm -f $(tests_util_containers_OBJS) Test $(tests_util_containers_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_containers_OBJS:.o=.d)

Test: Test.o $(util_LIB) Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(tests_util_containers_OBJS:.o=.d)

