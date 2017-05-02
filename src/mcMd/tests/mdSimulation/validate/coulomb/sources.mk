mcMd_tests_mdSimulation_validate_coulomb_= \
   mcMd/tests/mdSimulation/validate/coulomb/EwaldTest.cc

mcMd_tests_mdSimulation_validate_coulomb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_mdSimulation_validate_coulomb_))
mcMd_tests_mdSimulation_validate_coulomb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_mdSimulation_validate_coulomb_:.cc=.o))

