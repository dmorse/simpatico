SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(util_tests_accumulators_stochastic_autocorrstage_OBJS) AutoCorrStageTest

clean:
	rm -f $(util_tests_accumulators_stochastic_autocorrstage_OBJS) AutoCorrStageTest 
	rm -f $(util_tests_accumulators_stochastic_autocorrstage_OBJS:.o=.d)

clean-deps:
	rm -f $(util_tests_accumulators_stochastic_autocorrstage_OBJS:.o=.d)

AutoCorrStageTest: AutoCorrStageTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AutoCorrStageTest AutoCorrStageTest.o \
	       $(util_LIB)

-include $(util_tests_accumulators_stochastic_autocorrstage_OBJS:.o=.d)

