SRC_DIR_REL =../../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
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

