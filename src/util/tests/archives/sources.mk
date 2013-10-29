ifdef UTIL_MPI
util_tests_archives_=util/tests/archives/MpiTest.cpp
else
util_tests_archives_=util/tests/archives/Test.cpp
endif

util_tests_archives_SRCS=\
     $(addprefix $(SRC_DIR)/, $(util_tests_archives_))
util_tests_archives_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(util_tests_archives_:.cpp=.o))

