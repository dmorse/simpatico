ddMd_tests_sp_configIos_=ddMd/tests/sp/configIos/Test.cc

ddMd_tests_sp_configIos_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_sp_configIos_))
ddMd_tests_sp_configIos_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_sp_configIos_:.cc=.o))

