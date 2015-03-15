mcMd_tests_species_=mcMd/tests/species/Test.cc

mcMd_tests_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(mcMd_tests_species_))
mcMd_tests_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(mcMd_tests_species_:.cc=.o))

