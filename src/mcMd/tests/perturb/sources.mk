mcMd_tests_perturb_=mcMd/tests/perturb/Test.cc 

mcMd_tests_perturb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_perturb_))
mcMd_tests_perturb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_perturb_:.cc=.o))
mcMd_tests_perturb_EXES=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_perturb_:.cc=))

