SRC_DIR_REL =../../../../../src

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk


include sources.mk

all: $(tests_util_accumulators_stochastic_average_OBJS) AverageTest

clean:
	rm -f $(tests_util_accumulators_stochastic_average_OBJS) AverageTest $(tests_util_accumulators_stochastic_average_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_stochastic_average_OBJS:.o=.d)

AverageTest: AverageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageTest AverageTest.o \
	       $(util_LIB)

-include $(tests_util_accumulators_stochastic_average_OBJS:.o=.d)

