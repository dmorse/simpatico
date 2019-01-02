BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/mdPp/tests/storage/sources.mk

TEST=$(SRC_DIR)/mdPp/tests/storage/Test

all: $(mdPp_tests_storage_OBJS)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mdPp_tests_storage_OBJS) 
	rm -f $(mdPp_tests_storage_OBJS:.o=.d)
	rm -f $(mdPp_tests_storage_OBJS:.o=)

-include $(mdPp_tests_storage_OBJS:.o=.d)
-include $(mdPp_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

