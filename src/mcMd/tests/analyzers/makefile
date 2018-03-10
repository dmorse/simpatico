BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/analyzers/sources.mk

TEST=$(SRC_DIR)/mcMd/tests/analyzers/Test

.Secondary: $(mcMd_tests_analyzers_EXES) 

all: $(mcMd_tests_analyzers_EXES) 

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_analyzers_EXES)
	rm -f $(mcMd_tests_analyzers_OBJS)
	rm -f $(mcMd_tests_analyzers_OBJS:.o=.d)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f ar.txt binary

-include $(mcMd_tests_analyzers_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

