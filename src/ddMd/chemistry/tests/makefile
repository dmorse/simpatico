SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include sources.mk

all: $(ddMd_chemistry_tests_OBJS) AtomArrayTest GroupTest

clean:
	rm -f $(ddMd_chemistry_tests_OBJS) AtomArrayTest GroupTest $(ddMd_chemistry_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_chemistry_tests_OBJS:.o=.d)

AtomArrayTest: AtomArrayTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o AtomArrayTest AtomArrayTest.o \
	       $(ddMd_LIB) $(util_LIB)

GroupTest: GroupTest.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o GroupTest GroupTest.o \
	       $(ddMd_LIB) $(util_LIB)

-include $(ddMd_chemistry_tests_OBJS:.o=.d)

