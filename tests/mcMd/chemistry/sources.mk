tests_mcMd_chemistry_SRCS=$(TESTS_DIR)/mcMd/chemistry/MpiTest.cc \
    $(TESTS_DIR)/mcMd/chemistry/Test.cc 

tests_mcMd_chemistry_OBJS=$(tests_mcMd_chemistry_SRCS:.cc=.o)

