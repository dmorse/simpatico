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

ifdef UTIL_MPI
all: $(tests_mcMd_chemistry_OBJS) MpiTest 
else
all: $(tests_mcMd_chemistry_OBJS) Test
endif

clean:
	rm -f $(tests_mcMd_chemistry_OBJS) MpiTest Test $(tests_mcMd_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_mcMd_chemistry_OBJS:.o=.d)

MpiTest: MpiTest.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(mcMd_LIB) $(util_LIB)

Test: Test.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include $(tests_mcMd_chemistry_OBJS:.o=.d)
-include Test.d

