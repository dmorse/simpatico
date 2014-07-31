ddMd_tests_sp_neighbor_= ddMd/tests/sp/neighbor/Test.cc

ddMd_tests_sp_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_neighbor_))
ddMd_tests_sp_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_neighbor_:.cc=.o))

