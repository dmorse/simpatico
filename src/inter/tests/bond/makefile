OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/bond/sources.mk

TEST=inter/tests/bond/Test

all: $(inter_tests_bond_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(inter_tests_bond_OBJS) $(inter_tests_bond_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 

-include $(inter_tests_bond_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

