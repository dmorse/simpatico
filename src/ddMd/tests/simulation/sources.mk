ddMd_tests_simulation_=ddMd/tests/simulation/Test.cpp

ddMd_tests_simulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_simulation_))
ddMd_tests_simulation_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_simulation_:.cpp=.o))

