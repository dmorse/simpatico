mcMd_tests_neighbor_=mcMd/tests/neighbor/Test.cc

mcMd_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_neighbor_))
mcMd_tests_neighbor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_neighbor_:.cc=.o))

