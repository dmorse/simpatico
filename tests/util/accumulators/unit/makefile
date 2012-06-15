SRC_DIR_REL =../../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(tests_util_accumulators_unit_OBJS) Test

clean:
	rm -f Test $(tests_util_accumulators_unit_OBJS)
	rm -f $(tests_util_accumulators_unit_OBJS)
	rm -f $(tests_util_accumulators_unit_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_unit_OBJS:.o=.d)


Test: $(util_LIB) $(tests_util_accumulators_unit_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test \
        $(TESTS_DIR)/util/accumulators/unit/Test.o \
	       $(util_LIB)

-include $(tests_util_accumulators_unit_OBJS:.o=.d)

