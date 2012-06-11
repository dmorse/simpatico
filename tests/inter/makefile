SRC_DIR_REL =../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/inter/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(tests_inter_OBJS) Test

clean:
	rm -f $(tests_inter_OBJS) Test $(tests_inter_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_inter_OBJS:.o=.d)

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(tests_inter_OBJS:.o=.d)

