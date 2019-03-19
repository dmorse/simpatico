BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/analysis/sources.mk

TEST=$(SRC_DIR)/simp/tests/analysis/Test

.Secondary: $(simp_tests_analysis_OBJS)

all: $(simp_tests_analysis_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(simp_tests_analysis_OBJS) 
	rm -f $(simp_tests_analysis_EXES)
	rm -f $(simp_tests_analysis_OBJS:.o=.d)
	$(MAKE) clean-outputs

clean-outputs:
	@rm -f tmp/*

-include $(simp_tests_analysis_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

