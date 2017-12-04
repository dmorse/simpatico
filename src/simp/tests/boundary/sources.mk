simp_tests_boundary_=simp/tests/boundary/Test.cc

simp_tests_boundary_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_boundary_))
simp_tests_boundary_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_boundary_:.cc=.o))

