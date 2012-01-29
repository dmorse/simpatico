SRC_DIR_REL =../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/test/sources.mk
include sources.mk

all: $(util_param_tests_mpi_OBJS) MpiTest

clean:
	rm -f $(util_param_tests_mpi_OBJS) MpiTest $(util_param_tests_mpi_OBJS:.o=.d)

clean-deps:
	rm -f $(util_param_tests_mpi_OBJS:.o=.d)

MpiTest: MpiTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(util_LIB)

-include $(util_param_tests_mpi_OBJS:.o=.d)

