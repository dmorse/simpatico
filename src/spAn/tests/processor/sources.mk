spAn_tests_processor_=spAn/tests/processor/Test.cc

spAn_tests_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_processor_))
spAn_tests_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_processor_:.cc=.o))

