OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/potentials/sources.mk

all: $(ddMd_tests_potentials_OBJS) Test

clean:
	rm -f $(ddMd_tests_potentials_OBJS) Test $(ddMd_tests_potentials_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_tests_potentials_OBJS:.o=.d)

Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_potentials_OBJS:.o=.d)

