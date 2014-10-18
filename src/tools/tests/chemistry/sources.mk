tools_tests_chemistry_= tools/tests/chemistry/Test.cc

tools_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_chemistry_))
tools_tests_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_chemistry_:.cc=.o))

