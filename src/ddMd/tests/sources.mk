ddMd_tests_sp_=ddMd/tests/Test.cc

ddMd_tests_sp_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_))
ddMd_tests_sp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_:.cc=.o))

