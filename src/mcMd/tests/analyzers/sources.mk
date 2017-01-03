mcMd_tests_analyzers_=mcMd/tests/analyzers/Test.cc 

mcMd_tests_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_analyzers_))
mcMd_tests_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_analyzers_:.cc=.o))

