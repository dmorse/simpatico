mcMd_tests_simulation_=mcMd/tests/simulation/Test.cc

mcMd_tests_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_simulation_))
mcMd_tests_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_simulation_:.cc=.o))

