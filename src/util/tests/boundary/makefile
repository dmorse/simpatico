OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/boundary/sources.mk

TEST=util/tests/boundary/Test

all: $(util_tests_boundary_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(util_tests_boundary_OBJS) $(OBJ_DIR)/$(TEST) $(util_tests_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(util_tests_boundary_OBJS:.o=.d)

-include $(util_tests_boundary_OBJS:.o=.d)

