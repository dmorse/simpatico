spAn_tests_storage_=spAn/tests/storage/Test.cc

spAn_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(spAn_tests_storage_))
spAn_tests_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(spAn_tests_storage_:.cc=.o))

