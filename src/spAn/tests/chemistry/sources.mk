spAn_tests_chemistry_= spAn/tests/chemistry/Test.cc

spAn_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_chemistry_))
spAn_tests_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_chemistry_:.cc=.o))

