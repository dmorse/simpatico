SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include sources.mk

all: $(util_mpi_tests_OBJS) Test

clean:
	rm -f $(util_mpi_tests_OBJS) Test $(util_mpi_tests_OBJS:.o=.d)

clean-deps:
	rm -f $(util_mpi_tests_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_mpi_tests_OBJS:.o=.d)

