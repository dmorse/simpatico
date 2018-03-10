BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/tests/mdIntegrators/sources.mk
include $(SRC_DIR)/mcMd/patterns.mk

TEST=$(BLD_DIR)/mcMd/tests/mdIntegrators/Test

.Secondary: $(mcMd_tests_mdIntegrators_OBJS)

all: $(mcMd_tests_mdIntegrators_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_mdIntegrators_OBJS) 
	rm -f $(mcMd_tests_mdIntegrators_OBJS:.o=.d)
	rm -f $(mcMd_tests_mdIntegrators_OBJS:.o=)

-include $(mcMd_tests_mdIntegrators_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

