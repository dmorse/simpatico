BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/ddMd/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/chemistry/sources.mk

all: $(ddMd_tests_chemistry_OBJS) AtomTest GroupTest

clean:
	rm -f $(ddMd_tests_chemistry_OBJS) $(ddMd_tests_chemistry_OBJS:.o=.d)
	rm -f AtomTest GroupTest

clean-deps:
	rm -f $(ddMd_tests_chemistry_OBJS:.o=.d)

AtomTest: AtomTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AtomTest AtomTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

GroupTest: GroupTest.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o GroupTest GroupTest.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_chemistry_OBJS:.o=.d)

