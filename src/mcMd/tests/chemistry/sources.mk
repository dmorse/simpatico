ifdef UTIL_MPI
mcMd_tests_chemistry_=mcMd/tests/chemistry/MpiTest.cc 
else
mcMd_tests_chemistry_=mcMd/tests/chemistry/Test.cc 
endif

mcMd_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_chemistry_))
mcMd_tests_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_chemistry_:.cc=.o))

