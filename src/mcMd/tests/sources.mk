mcMd_tests_=mcMd/tests/Test.cc

mcMd_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_))
mcMd_tests_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_:.cc=.o))

