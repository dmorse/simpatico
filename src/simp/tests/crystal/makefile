BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/crystal/sources.mk

TEST=simp/tests/crystal/Test

all: $(simp_tests_crystal_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(simp_tests_crystal_OBJS) $(BLD_DIR)/$(TEST) $(simp_tests_crystal_OBJS:.o=.d)

-include $(simp_tests_crystal_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)

