mdCf_tests_storage_=mdCf/tests/storage/Test.cc

mdCf_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdCf_tests_storage_))
mdCf_tests_storage_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdCf_tests_storage_:.cc=.o))

