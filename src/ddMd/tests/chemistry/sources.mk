tests_ddMd_chemistry_SRCS=$(TESTS_DIR)/ddMd/chemistry/AtomTest.cc \
    $(TESTS_DIR)/ddMd/chemistry/GroupTest.cc 

tests_ddMd_chemistry_OBJS=$(tests_ddMd_chemistry_SRCS:.cc=.o)

ddMd_tests_TMP_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_TMP_))
ddMd_tests_TMP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_TMP_:.cpp=.o))

