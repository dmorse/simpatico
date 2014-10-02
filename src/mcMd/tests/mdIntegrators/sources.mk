mcMd_tests_mdIntegrators_=mcMd/tests/mdIntegrators/Test.cc

mcMd_tests_mdIntegrators_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_mdIntegrators_))
mcMd_tests_mdIntegrators_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_mdIntegrators_:.cc=.o))

