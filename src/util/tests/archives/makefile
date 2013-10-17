SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk

include sources.mk

ifdef UTIL_MPI

all: $(tests_util_archives_OBJS) MpiTest

clean:
	rm -f $(tests_util_archives_OBJS) MpiTest $(tests_util_archives_OBJS:.o=.d)
	rm -f text binary dummy

MpiTest: $(util_LIB) $(tests_util_archives_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(util_LIB)
else

all: $(tests_util_archives_OBJS) Test

clean:
	rm -f $(tests_util_archives_OBJS) Test $(tests_util_archives_OBJS:.o=.d)
	rm -f text binary dummy

Test: $(util_LIB) $(tests_util_archives_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

endif

clean-deps:
	rm -f $(tests_util_archives_OBJS:.o=.d)

-include $(tests_util_archives_OBJS:.o=.d)

