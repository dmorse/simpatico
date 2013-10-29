OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/dihedral/sources.mk

all: $(inter_tests_dihedral_OBJS) Test

clean:
	rm -f $(inter_tests_dihedral_OBJS) Test $(inter_tests_dihedral_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_tests_dihedral_OBJS:.o=.d)

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(inter_tests_dihedral_OBJS:.o=.d)

