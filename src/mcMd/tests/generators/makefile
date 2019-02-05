BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/tests/generators/sources.mk

TEST=$(BLD_DIR)/mcMd/tests/generators/Test
.Secondary: $(mcMd_tests_generators_OBJS)

all: $(mcMd_tests_generators_EXES)

run: $(TEST)
	$(TEST)

clean:
	rm -f $(mcMd_tests_generators_OBJS) 
	rm -f $(mcMd_tests_generators_OBJS:.o=.d)
	rm -f $(mcMd_tests_generators_OBJS:.o=)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f out/*
	rm -f *.rst

-include $(mcMd_tests_generators_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

