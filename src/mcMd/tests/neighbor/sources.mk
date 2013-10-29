mcMd_tests_neighbor_=mcMd/tests/neighbor/Test.cpp

mcMd_tests_neighbor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_neighbor_))
mcMd_tests_neighbor_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_neighbor_:.cpp=.o))

