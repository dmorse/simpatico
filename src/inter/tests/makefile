SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include sources.mk

all: $(inter_tests_OBJS) Test

clean:
	rm -f $(inter_tests_OBJS) Test $(inter_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_tests_OBJS:.o=.d)

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(inter_tests_OBJS:.o=.d)

