OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/format/sources.mk

TEST=util/tests/format/TEST

all: $(util_tests_format_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(util_tests_format_OBJS) $(util_tests_format_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 

-include $(util_tests_format_OBJS:.o=.d)

