SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(util_accumulators_tests_stochastic_averageStage_OBJS) AverageStageTest

clean:
	rm -f $(util_accumulators_tests_stochastic_averageStage_OBJS) AverageStageTest $(util_accumulators_tests_stochastic_averageStage_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_tests_stochastic_averageStage_OBJS:.o=.d)

AverageStageTest: AverageStageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AverageStageTest AverageStageTest.o \
	       $(util_LIB)

-include $(util_accumulators_tests_stochastic_averageStage_OBJS:.o=.d)

