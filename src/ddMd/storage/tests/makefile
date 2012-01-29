SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/ddMd/defines.mk
include $(SRC_DIR_REL)/ddMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/ddMd/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(ddMd_storage_tests_OBJS) Test

clean:
	rm -f $(ddMd_storage_tests_OBJS) Test $(ddMd_storage_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_storage_tests_OBJS:.o=.d)

Test: Test.o $(ddMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(util_LIB)

-include $(ddMd_storage_tests_OBJS:.o=.d)

