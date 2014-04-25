mdPp_tests_processor_=mdPp/tests/processor/Test.cc

mdPp_tests_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_processor_))
mdPp_tests_processor_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_tests_processor_:.cc=.o))

