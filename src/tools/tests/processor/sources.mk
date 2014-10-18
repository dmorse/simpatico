tools_tests_processor_=tools/tests/processor/Test.cc

tools_tests_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_processor_))
tools_tests_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_processor_:.cc=.o))

