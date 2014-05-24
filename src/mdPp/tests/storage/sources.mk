mdPp_tests_storage_=mdPp/tests/storage/Test.cc

mdPp_tests_storage_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mdPp_tests_storage_))
mdPp_tests_storage_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mdPp_tests_storage_:.cc=.o))

