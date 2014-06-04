ddMd_tests_sp_chemistry_= mdCf/tests/chemistry/Test.cc

ddMd_tests_sp_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_chemistry_))
ddMd_tests_sp_chemistry_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_tests_sp_chemistry_:.cc=.o))

