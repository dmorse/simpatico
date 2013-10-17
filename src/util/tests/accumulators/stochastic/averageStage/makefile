SRC_DIR_REL =../../../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk


include sources.mk

all: $(tests_util_accumulators_stochastic_averageStage_OBJS) AverageStageTest

clean:
	rm -f $(tests_util_accumulators_stochastic_averageStage_OBJS) AverageStageTest $(tests_util_accumulators_stochastic_averageStage_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_stochastic_averageStage_OBJS:.o=.d)

AverageStageTest: AverageStageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageStageTest AverageStageTest.o \
	       $(util_LIB)

-include $(tests_util_accumulators_stochastic_averageStage_OBJS:.o=.d)

