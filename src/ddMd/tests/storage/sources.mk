tests_ddMd_storage_SRCS=$(TESTS_DIR)/ddMd/storage/Test.cc 

tests_ddMd_storage_OBJS=$(tests_ddMd_storage_SRCS:.cc=.o)

ddMd_tests_TMP_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_TMP_))
ddMd_tests_TMP_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_TMP_:.cpp=.o))

