SRC_DIR_REL =../../../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(tests_util_accumulators_stochastic_autocorr_OBJS) AutoCorrTest

clean:
	rm -f $(tests_util_accumulators_stochastic_autocorr_OBJS) AutoCorrTest $(tests_util_accumulators_stochastic_autocorr_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_stochastic_autocorr_OBJS:.o=.d)

AutoCorrTest: AutoCorrTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AutoCorrTest AutoCorrTest.o \
	       $(util_LIB)

-include $(tests_util_accumulators_stochastic_autocorr_OBJS:.o=.d)

