SRC_DIR_REL =../../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk


include sources.mk

all: $(tests_util_accumulators_stochastic_OBJS)
	cd autocorr; make all
	cd average;  make all
	cd averageStage; make all
	cd meanSqDisp; make all
	cd random; make all

clean:
	cd autocorr; make clean
	cd average;  make clean
	cd averageStage; make clean
	cd meanSqDisp; make clean
	cd random; make clean
	rm -f $(tests_util_accumulators_stochastic_OBJS) $(tests_util_accumulators_stochastic_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_accumulators_stochastic_OBJS:.o=.d)

-include $(tests_util_accumulators_stochastic_OBJS:.o=.d)

