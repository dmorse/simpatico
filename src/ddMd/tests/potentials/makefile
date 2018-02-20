BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/potentials/sources.mk

TEST= $(SRC_DIR)/ddMd/tests/potentials/Test
.Secondary: $(ddMd_tests_potentials_OBJS)

all: $(ddMd_tests_potentials_EXES)

run: $(TEST)
	$(MPIRUN) 6 $(TEST)

clean:
	rm -f $(ddMd_tests_potentials_EXES)
	rm -f $(ddMd_tests_potentials_OBJS) 
	rm -f $(ddMd_tests_potentials_OBJS:.o=.d)

-include $(ddMd_tests_potentials_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

