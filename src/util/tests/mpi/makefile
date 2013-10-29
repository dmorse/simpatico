OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/TMP/sources.mk

all: $(util_tests_TMP_OBJS) Test

clean:
	rm -f $(util_tests_TMP_OBJS) Test $(util_tests_TMP_OBJS:.o=.d)

clean-deps:
	rm -f $(util_tests_TMP_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_tests_TMP_OBJS:.o=.d)

