SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/ddMd/defines.mk
include $(TESTS_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/sources.mk


include sources.mk

all: $(tests_ddMd_chemistry_OBJS) AtomArrayTest GroupTest

clean:
	rm -f $(tests_ddMd_chemistry_OBJS) AtomArrayTest GroupTest $(tests_ddMd_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_ddMd_chemistry_OBJS:.o=.d)

AtomArrayTest: AtomArrayTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AtomArrayTest AtomArrayTest.o \
	       $(ddMd_LIB) $(util_LIB)

GroupTest: GroupTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o GroupTest GroupTest.o \
	       $(ddMd_LIB) $(util_LIB)

-include $(tests_ddMd_chemistry_OBJS:.o=.d)

