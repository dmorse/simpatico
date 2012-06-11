SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(tests_mcMd_links_OBJS) Test

clean:
	rm -f $(tests_mcMd_links_OBJS) Test $(tests_mcMd_links_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_mcMd_links_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include $(tests_mcMd_links_OBJS:.o=.d)

