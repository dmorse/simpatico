simp_tests_=simp/tests/Test.cc

simp_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_))
simp_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_:.cc=.o))

