ifdef UTIL_MPI
util_archives_tests_SRCS=$(SRC_DIR)/util/archives/tests/MpiTest.cc
else
util_archives_tests_SRCS=$(SRC_DIR)/util/archives/tests/Test.cc 
endif

util_archives_tests_OBJS=$(util_archives_tests_SRCS:.cc=.o)

