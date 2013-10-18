ddMd_tests_=ddMd/tests/Test.cpp

ddMd_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_))
ddMd_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_:.cpp=.o))

