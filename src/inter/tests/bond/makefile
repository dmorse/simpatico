BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/bond/sources.mk

all: $(inter_tests_bond_OBJS) Test

clean:
	rm -f $(inter_tests_bond_OBJS) Test $(inter_tests_bond_OBJS:.o=.d)

clean-deps:
	rm -f $(inter_tests_bond_OBJS:.o=.d)

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(inter_tests_bond_OBJS:.o=.d)

