BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/simulation/sources.mk

TEST= $(SRC_DIR)/mcMd/tests/simulation/Test

.Secondary: $(mcMd_tests_simulation_OBJS)

all: $(mcMd_tests_simulation_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_simulation_EXES)
	rm -f $(mcMd_tests_simulation_OBJS) 
	rm -f $(mcMd_tests_simulation_OBJS:.o=.d)

-include $(mcMd_tests_simulation_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

