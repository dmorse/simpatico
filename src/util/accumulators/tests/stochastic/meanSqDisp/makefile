SRC_DIR_REL =../../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(util_accumulators_tests_stochastic_meanSqDisp_OBJS) MeanSqDispTest

clean:
	rm -f $(util_accumulators_tests_stochastic_meanSqDisp_OBJS) MeanSqDispTest $(util_accumulators_tests_stochastic_meanSqDisp_OBJS:.o=.d)

clean-deps:
	rm -f $(util_accumulators_tests_stochastic_meanSqDisp_OBJS:.o=.d)

MeanSqDispTest: MeanSqDispTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MeanSqDispTest MeanSqDispTest.o \
	       $(util_LIB)

-include $(util_accumulators_tests_stochastic_meanSqDisp_OBJS:.o=.d)

