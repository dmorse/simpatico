mcMd_tests_potentials_coulomb_= \
   mcMd/tests/potentials/coulomb/EwaldTest.cc \
   mcMd/tests/potentials/coulomb/SpmeTest.cc

mcMd_tests_potentials_coulomb_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_potentials_coulomb_))
mcMd_tests_potentials_coulomb_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_potentials_coulomb_:.cc=.o))

