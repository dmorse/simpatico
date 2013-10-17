SRC_DIR_REL =../../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(tests_util_random_OBJS) Test

clean:
	rm -f $(tests_util_random_OBJS) Test $(tests_util_random_OBJS:.o=.d) binary

clean-deps:
	rm -f $(tests_util_random_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(tests_util_random_OBJS:.o=.d)

