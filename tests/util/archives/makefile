
SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

ifdef UTIL_MPI

all: $(tests_util_archives_OBJS) MpiTest

clean:
	rm -f $(tests_util_archives_OBJS) MpiTest $(tests_util_archives_OBJS:.o=.d)

MpiTest: $(util_LIB) $(tests_util_archives_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(util_LIB)
else

all: $(tests_util_archives_OBJS) Test

clean:
	rm -f $(tests_util_archives_OBJS) Test $(tests_util_archives_OBJS:.o=.d)

Test: $(util_LIB) $(tests_util_archives_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

endif

clean-deps:
	rm -f $(tests_util_archives_OBJS:.o=.d)

-include $(tests_util_archives_OBJS:.o=.d)

