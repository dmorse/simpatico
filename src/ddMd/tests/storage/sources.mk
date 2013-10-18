ddMd_tests_storage_=ddMd/tests/storage/Test.cpp

ddMd_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_storage_))
ddMd_tests_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_storage_:.cpp=.o))

