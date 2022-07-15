simp_tests_crystal_=simp/tests/crystal/Test.cc

simp_tests_crystal_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_crystal_))
simp_tests_crystal_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_crystal_:.cc=.o))

