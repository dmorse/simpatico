ddMd_tests_=ddMd/tests/Test.cc

ddMd_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_))
ddMd_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_:.cc=.o))
ddMd_tests_EXES=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_:.cc=))

