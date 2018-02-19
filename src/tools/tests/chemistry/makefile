BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/tools/tests/chemistry/sources.mk

TEST= $(SRC_DIR)/tools/tests/chemistry/Test

all: $(tools_tests_chemistry_OBJS)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(tools_tests_chemistry_OBJS) 
	rm -f $(tools_tests_chemistry_OBJS:.o=.d)
	rm -f $(tools_tests_chemistry_OBJS:.o=)

-include $(tools_tests_chemistry_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

