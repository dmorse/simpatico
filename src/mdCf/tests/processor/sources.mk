mdCf_tests_processor_=mdCf/tests/processor/Test.cc

mdCf_tests_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_tests_processor_))
mdCf_tests_processor_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_tests_processor_:.cc=.o))

