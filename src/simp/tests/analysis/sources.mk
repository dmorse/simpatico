simp_tests_analysis_=simp/tests/analysis/Test.cc

simp_tests_analysis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_analysis_))
simp_tests_analysis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_analysis_:.cc=.o))
simp_tests_analysis_EXES=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_analysis_:.cc=))

