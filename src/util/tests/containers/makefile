BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/containers/sources.mk

all: $(util_tests_containers_OBJS) Test

clean:
	rm -f $(util_tests_containers_OBJS) Test $(util_tests_containers_OBJS:.o=.d)

clean-deps:
	rm -f $(util_tests_containers_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_tests_containers_OBJS:.o=.d)

