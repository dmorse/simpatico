BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/mcSimulation/sources.mk

TEST=$(BLD_DIR)/mcMd/tests/mcSimulation/Test

.Secondary: $(mcMd_tests_mcSimulation_OBJS)

all: $(mcMd_tests_mcSimulation_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_mcSimulation_OBJS) 
	rm -f $(mcMd_tests_mcSimulation_OBJS:.o=.d)
	rm -f $(mcMd_tests_mcSimulation_EXES)
	$(MAKE) clean-outputs

clean-outputs:
	@rm -f tmp/*

-include $(mcMd_tests_mcSimulation_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

