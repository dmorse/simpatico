mcMd_tests_generators_=mcMd/tests/generators/Test.cc

mcMd_tests_generators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_generators_))
mcMd_tests_generators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_generators_:.cc=.o))

