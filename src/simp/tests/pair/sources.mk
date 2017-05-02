simp_tests_pair_=simp/tests/pair/Test.cc

simp_tests_pair_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_pair_))
simp_tests_pair_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_pair_:.cc=.o))

