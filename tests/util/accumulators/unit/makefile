SRC_DIR_REL =../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

all: $(util_accumulators_tests_unit_OBJS) Test

clean:
	rm -f Test $(util_accumulators_tests_unit_OBJS) \
              $(util_accumulators_tests_unit_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_tests_unit_OBJS:.o=.d)


Test: $(util_LIB) $(util_accumulators_tests_unit_OBJS)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_accumulators_tests_unit_OBJS:.o=.d)

