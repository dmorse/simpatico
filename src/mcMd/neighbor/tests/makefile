SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include sources.mk

all: $(mcMd_neighbor_tests_OBJS) Test

clean:
	rm -f $(mcMd_neighbor_tests_OBJS) Test $(mcMd_neighbor_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_neighbor_tests_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include $(mcMd_neighbor_tests_OBJS:.o=.d)

