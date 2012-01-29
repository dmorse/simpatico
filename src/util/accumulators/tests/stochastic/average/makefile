SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(util_accumulators_tests_stochastic_average_OBJS) AverageTest

clean:
	rm -f $(util_accumulators_tests_stochastic_average_OBJS) AverageTest $(util_accumulators_tests_stochastic_average_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_tests_stochastic_average_OBJS:.o=.d)

AverageTest: AverageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageTest AverageTest.o \
	       $(util_LIB)

-include $(util_accumulators_tests_stochastic_average_OBJS:.o=.d)

