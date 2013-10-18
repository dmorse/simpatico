ddMd_tests_neighbor_=ddMd/tests/neighbor/Test.cpp

ddMd_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_neighbor_))
ddMd_tests_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_neighbor_:.cpp=.o))

