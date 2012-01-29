SRC_DIR_REL =../../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/util/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include sources.mk

all: $(util_param_tests_serial_OBJS) Test

clean:
	rm -f $(util_param_tests_serial_OBJS) Test $(util_param_tests_serial_OBJS:.o=.d)

clean-deps:
	rm -f $(util_param_tests_serial_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_param_tests_serial_OBJS:.o=.d)

