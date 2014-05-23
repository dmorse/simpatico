mdPp_tests_chemistry_=mdPp/tests/chemistry/Test.cc

mdPp_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_chemistry_))
mdPp_tests_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_tests_chemistry_:.cc=.o))

