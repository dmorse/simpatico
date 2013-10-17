tests_ddMd_simulation_SRCS=$(TESTS_DIR)/ddMd/simulation/Test.cc 

tests_ddMd_simulation_OBJS=$(tests_ddMd_simulation_SRCS:.cc=.o)

ddMd_tests_TMP_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_TMP_))
ddMd_tests_TMP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_TMP_:.cpp=.o))

