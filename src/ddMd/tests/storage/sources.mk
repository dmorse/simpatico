ddMd_tests_storage_=ddMd/tests/storage/Test.cc

ddMd_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_storage_))
ddMd_tests_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_storage_:.cc=.o))

