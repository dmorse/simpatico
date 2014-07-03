ddMd_tests_sp_storage_=ddMd/tests/sp/storage/Test.cc

ddMd_tests_sp_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_storage_))
ddMd_tests_sp_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_storage_:.cc=.o))

