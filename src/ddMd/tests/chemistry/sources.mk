ddMd_tests_chemistry_=\
    ddMd/tests/chemistry/AtomTest.cc \
    ddMd/tests/chemistry/GroupTest.cc 

ddMd_tests_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_chemistry_))
ddMd_tests_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_tests_chemistry_:.cc=.o))

