SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

all: $(util_containers_tests_OBJS) Test

clean:
	rm -f $(util_containers_tests_OBJS) Test $(util_containers_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(util_containers_tests_OBJS:.o=.d)

Test: Test.o $(util_LIB) Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_containers_tests_OBJS:.o=.d)

