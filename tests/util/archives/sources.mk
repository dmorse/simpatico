ifdef UTIL_MPI
tests_util_archives_SRCS=$(TESTS_DIR)/util/archives/MpiTest.cc
else
tests_util_archives_SRCS=$(TESTS_DIR)/util/archives/Test.cc 
endif

tests_util_archives_OBJS=$(tests_util_archives_SRCS:.cc=.o)

