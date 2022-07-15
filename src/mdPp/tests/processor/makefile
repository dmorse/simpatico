BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/mdPp/tests/processor/sources.mk

TEST= $(SRC_DIR)/mdPp/tests/processor/Test

all: $(mdPp_tests_processor_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mdPp_tests_processor_OBJS) 
	rm -f $(mdPp_tests_processor_OBJS:.o=.d)
	rm -f $(mdPp_tests_processor_OBJS:.o=)
	rm -f out/*

-include $(mdPp_tests_processor_OBJS:.o=.d)
-include $(mdPp_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

