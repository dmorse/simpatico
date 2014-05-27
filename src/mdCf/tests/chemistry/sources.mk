mdCf_tests_chemistry_= mdCf/tests/chemistry/Test.cc

mdCf_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_tests_chemistry_))
mdCf_tests_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_tests_chemistry_:.cc=.o))

