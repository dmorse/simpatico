ifdef UTIL_MPI
mcMd_tests_chemistry_=mcMd/tests/chemistry/MpiTest.cpp 
else
mcMd_tests_chemistry_=mcMd/tests/chemistry/Test.cpp 
endif

mcMd_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_chemistry_))
mcMd_tests_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_chemistry_:.cpp=.o))

