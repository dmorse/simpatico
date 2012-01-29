SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(util_accumulators_tests_stochastic_autocorr_OBJS) AutoCorrTest

clean:
	rm -f $(util_accumulators_tests_stochastic_autocorr_OBJS) AutoCorrTest $(util_accumulators_tests_stochastic_autocorr_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_tests_stochastic_autocorr_OBJS:.o=.d)

AutoCorrTest: AutoCorrTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AutoCorrTest AutoCorrTest.o \
	       $(util_LIB)

-include $(util_accumulators_tests_stochastic_autocorr_OBJS:.o=.d)

