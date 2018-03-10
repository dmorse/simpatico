BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/tests/mdSimulation/sources.mk
include $(SRC_DIR)/mcMd/patterns.mk

TEST= $(BLD_DIR)/mcMd/tests/mdSimulation/Test

.Secondary: $(mcMd_tests_mdSimulation_OBJS)

all: $(mcMd_tests_mdSimulation_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_mdSimulation_OBJS) 
	rm -f $(mcMd_tests_mdSimulation_OBJS:.o=.d)
	rm -f $(mcMd_tests_mdSimulation_EXES)
	$(MAKE) clean-outputs

clean-outputs:
	@rm -f tmp/*

-include $(mcMd_tests_mdSimulation_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

