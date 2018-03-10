BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/simulation/sources.mk

TEST= $(SRC_DIR)/ddMd/tests/simulation/Test

.Secondary: $(ddMd_tests_simulation_OBJS) 

all: $(ddMd_tests_simulation_EXES) 

run: $(TEST)
	$(MPIRUN) 6 $(TEST)

clean:
	rm -f $(ddMd_tests_simulation_EXES)
	rm -f $(ddMd_tests_simulation_OBJS) 
	rm -f $(ddMd_tests_simulation_OBJS:.o=.d)

-include $(ddMd_tests_simulation_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

