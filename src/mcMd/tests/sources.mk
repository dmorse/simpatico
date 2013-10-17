mcMd_tests_=mcMd/tests/Test.cpp

mcMd_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_))
mcMd_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_:.cpp=.o))

