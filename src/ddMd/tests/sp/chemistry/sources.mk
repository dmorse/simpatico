ddMd_tests_sp_chemistry_= ddMd/tests/sp/chemistry/Test.cc

ddMd_tests_sp_chemistry_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_chemistry_))
ddMd_tests_sp_chemistry_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_chemistry_:.cc=.o))

