SRC_DIR_REL =../../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk


include sources.mk

all: $(tests_util_accumulators_unit_OBJS) Test

clean:
	-rm -f Test $(tests_util_accumulators_unit_OBJS)
	-rm -f $(tests_util_accumulators_unit_OBJS)
	-rm -f $(tests_util_accumulators_unit_OBJS:.o=.d)
	-rm -f binary

clean-deps:
	rm -f $(tests_util_accumulators_unit_OBJS:.o=.d)


Test: $(util_LIB) $(tests_util_accumulators_unit_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test \
        $(TESTS_DIR)/util/accumulators/unit/Test.o \
	       $(util_LIB)

-include $(tests_util_accumulators_unit_OBJS:.o=.d)

