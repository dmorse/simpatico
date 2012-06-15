include $(TESTS_DIR)/util/accumulators/stochastic/autocorr/sources.mk
include $(TESTS_DIR)/util/accumulators/stochastic/average/sources.mk
include $(TESTS_DIR)/util/accumulators/stochastic/averageStage/sources.mk
include $(TESTS_DIR)/util/accumulators/stochastic/meanSqDisp/sources.mk
include $(TESTS_DIR)/util/accumulators/stochastic/random/sources.mk

tests_util_accumulators_stochastic_SRCS=\
    $(tests_util_accumulators_stochastic_autocorr_SRCS) \
    $(tests_util_accumulators_stochastic_average_SRCS) \
    $(tests_util_accumulators_stochastic_averageStage_SRCS) \
    $(tests_util_accumulators_stochastic_meanSqDisp_SRCS) \
    $(tests_util_accumulators_stochastic_random_SRCS) 

tests_util_accumulators_stochastic_OBJS=$(tests_util_accumulators_stochastic_SRCS:.cc=.o)

