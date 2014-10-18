tools_tests_storage_=tools/tests/storage/Test.cc

tools_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(tools_tests_storage_))
tools_tests_storage_OBJS=\
     $(addprefix $(BLD_DIR)/, $(tools_tests_storage_:.cc=.o))

