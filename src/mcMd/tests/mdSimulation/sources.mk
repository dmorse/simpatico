mcMd_tests_mdSimulation_=mcMd/tests/mdSimulation/Test.cpp

mcMd_tests_mdSimulation_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_mdSimulation_))
mcMd_tests_mdSimulation_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(mcMd_tests_mdSimulation_:.cpp=.o))

