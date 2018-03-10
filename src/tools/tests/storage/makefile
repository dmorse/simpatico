BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/tools/tests/storage/sources.mk

TEST=$(SRC_DIR)/tools/tests/storage/Test

all: $(tools_tests_storage_OBJS)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(tools_tests_storage_OBJS) 
	rm -f $(tools_tests_storage_OBJS:.o=.d)
	rm -f $(tools_tests_storage_OBJS:.o=)

-include $(tools_tests_storage_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

