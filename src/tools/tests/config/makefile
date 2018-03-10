BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/tools/tests/config/sources.mk

TEST= $(SRC_DIR)/tools/tests/config/Test

all: $(tools_tests_config_OBJS)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(tools_tests_config_OBJS) 
	rm -f $(tools_tests_config_OBJS:.o=.d)
	rm -f $(tools_tests_config_OBJS:.o=)
	rm -f out/*

-include $(tools_tests_config_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

