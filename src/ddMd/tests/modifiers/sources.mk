ddMd_tests_modifiers_=ddMd/tests/modifiers/Test.cc

ddMd_tests_modifiers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_tests_modifiers_))
ddMd_tests_modifiers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_tests_modifiers_:.cc=.o))

