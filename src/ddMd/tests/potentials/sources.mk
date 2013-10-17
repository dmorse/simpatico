tests_ddMd_potentials_SRCS=$(TESTS_DIR)/ddMd/potentials/Test.cc 

tests_ddMd_potentials_OBJS=$(tests_ddMd_potentials_SRCS:.cc=.o)

ddMd_tests_TMP_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_TMP_))
ddMd_tests_TMP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_TMP_:.cpp=.o))

