mdPp_tests_processor_=mdPp/tests/processor/Test.cc

mdPp_tests_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_processor_))
mdPp_tests_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_processor_:.cc=.o))
mdPp_tests_processor_EXES=\
     $(addprefix $(BLD_DIR)/, $(mdPp_tests_processor_:.cc=))

