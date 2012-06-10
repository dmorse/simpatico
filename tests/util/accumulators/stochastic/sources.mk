include $(SRC_DIR)/util/accumulators/tests/stochastic/autocorr/sources.mk
include $(SRC_DIR)/util/accumulators/tests/stochastic/average/sources.mk
include $(SRC_DIR)/util/accumulators/tests/stochastic/averageStage/sources.mk
include $(SRC_DIR)/util/accumulators/tests/stochastic/meanSqDisp/sources.mk
include $(SRC_DIR)/util/accumulators/tests/stochastic/random/sources.mk

util_accumulators_tests_stochastic_SRCS=\
    $(util_accumulators_tests_stochastic_autocorr_SRCS) \
    $(util_accumulators_tests_stochastic_average_SRCS) \
    $(util_accumulators_tests_stochastic_averageStage_SRCS) \
    $(util_accumulators_tests_stochastic_meanSqDisp_SRCS) \
    $(util_accumulators_tests_stochastic_random_SRCS) 

util_accumulators_tests_stochastic_OBJS=$(util_accumulators_tests_stochastic_SRCS:.cc=.o)

