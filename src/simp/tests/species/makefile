BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/species/sources.mk

TEST=$(SRC_DIR)/simp/tests/species/Test

.Secondary: $(simp_tests_species_OBJS)

all: $(simp_tests_species_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(simp_tests_species_OBJS) 
	rm -f $(simp_tests_species_EXES)
	rm -f $(simp_tests_species_OBJS:.o=.d)
	$(MAKE) clean-outputs

clean-outputs:
	@rm -f tmp/*
	@rm -f binary
	@rm -f structure

-include $(simp_tests_species_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

