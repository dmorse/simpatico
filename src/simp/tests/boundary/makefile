BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/simp/tests/boundary/sources.mk

TEST=simp/tests/boundary/Test

all: $(simp_tests_boundary_OBJS) $(BLD_DIR)/$(TEST)

clean:
	rm -f $(simp_tests_boundary_OBJS) $(BLD_DIR)/$(TEST) $(simp_tests_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_tests_boundary_OBJS:.o=.d)

-include $(simp_tests_boundary_OBJS:.o=.d)

