SRC_DIR_REL =../../../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include sources.mk

all: $(tests_util_param_mpi_OBJS) MpiTest

clean:
	rm -f $(tests_util_param_mpi_OBJS) MpiTest $(tests_util_param_mpi_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_util_param_mpi_OBJS:.o=.d)

MpiTest: MpiTest.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o MpiTest MpiTest.o \
	       $(util_LIB)

-include $(tests_util_param_mpi_OBJS:.o=.d)

