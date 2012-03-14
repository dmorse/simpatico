SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk
include sources.mk

ifdef UTIL_MPI
all: $(mcMd_chemistry_tests_OBJS) MpiTest 
else
all: $(mcMd_chemistry_tests_OBJS) Test
endif

clean:
	rm -f $(mcMd_chemistry_tests_OBJS) MpiTest Test $(mcMd_chemistry_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_chemistry_tests_OBJS:.o=.d)

MpiTest: MpiTest.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(mcMd_LIB) $(util_LIB)

Test: Test.o $(mcMd_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(util_LIB)

-include $(mcMd_chemistry_tests_OBJS:.o=.d)
-include Test.d

