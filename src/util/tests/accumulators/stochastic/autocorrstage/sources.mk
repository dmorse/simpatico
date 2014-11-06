util_tests_accumulators_stochastic_autocorrstage_ = \
    util/tests/accumulators/stochastic/autocorrstage/AutoCorrStageTest.cc 

util_tests_accumulators_stochastic_autocorrstage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_accumulators_stochastic_autocorrstage_))
util_tests_accumulators_stochastic_autocorrstage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(util_tests_accumulators_stochastic_autocorrstage_:.cc=.o))

