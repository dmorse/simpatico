mcMd_tests_mcSimulation_=mcMd/tests/mcSimulation/Test.cc

mcMd_tests_mcSimulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_mcSimulation_))
mcMd_tests_mcSimulation_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_mcSimulation_:.cc=.o))

