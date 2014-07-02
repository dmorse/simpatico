ddMd_tests_sp_processor_=ddMd/tests/sp/processor/Test.cc

ddMd_tests_sp_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_processor_))
ddMd_tests_sp_processor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_processor_:.cc=.o))

