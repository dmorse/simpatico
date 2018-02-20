BLD_DIR_REL=../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/modifiers/sources.mk

TEST= $(SRC_DIR)/ddMd/tests/modifiers/Test

.Secondary: $(ddMd_tests_modifiers_OBJS)

all: $(ddMd_tests_modifiers_EXES)

run: $(TEST)
	$(MPIRUN) 6 $(TEST)

clean:
	rm -f $(ddMd_tests_modifiers_OBJS) 
	rm -f $(ddMd_tests_modifiers_OBJS:.o=.d)
	rm -f $(ddMd_tests_modifiers_OBJS:.o=)

clean-deps:
	rm -f $(ddMd_tests_modifiers_OBJS:.o=.d)

-include $(ddMd_tests_modifiers_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

