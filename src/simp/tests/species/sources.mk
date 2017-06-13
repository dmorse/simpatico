simp_tests_species_=simp/tests/species/Test.cc

simp_tests_species_SRCS=\
     $(addprefix $(SRC_DIR)/, $(simp_tests_species_))
simp_tests_species_OBJS=\
     $(addprefix $(BLD_DIR)/, $(simp_tests_species_:.cc=.o))

