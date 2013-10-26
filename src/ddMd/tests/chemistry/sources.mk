ddMd_tests_chemistry_=\
    ddMd/tests/chemistry/AtomTest.cpp \
    ddMd/tests/chemistry/GroupTest.cpp 

ddMd_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_chemistry_))
ddMd_tests_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_tests_chemistry_:.cpp=.o))

