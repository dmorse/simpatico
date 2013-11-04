OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/dihedral/sources.mk

TEST=inter/tests/dihedral/Test

all: $(inter_tests_dihedral_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(inter_tests_dihedral_OBJS) $(inter_tests_dihedral_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 

-include $(inter_tests_dihedral_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

