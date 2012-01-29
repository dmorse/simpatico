SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

ifdef UTIL_MPI

all: $(util_archives_tests_OBJS) MpiTest

clean:
	rm -f $(util_archives_tests_OBJS) MpiTest $(util_archives_tests_OBJS:.o=.d)

MpiTest: $(util_LIB) $(util_archives_tests_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(util_LIB)
else

all: $(util_archives_tests_OBJS) Test

clean:
	rm -f $(util_archives_tests_OBJS) Test $(util_archives_tests_OBJS:.o=.d)

Test: $(util_LIB) $(util_archives_tests_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

endif

clean-deps:
	rm -f $(util_archives_tests_OBJS:.o=.d)

-include $(util_archives_tests_OBJS:.o=.d)

