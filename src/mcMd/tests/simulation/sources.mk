mcMd_tests_simulation_=mcMd/tests/simulation/Test.cpp

mcMd_tests_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_simulation_))
mcMd_tests_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_simulation_:.cpp=.o))

