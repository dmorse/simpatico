SRC_DIR_REL =../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/ddMd/defines.mk
include $(TESTS_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

include sources.mk

all: $(tests_ddMd_potentials_OBJS) Test

clean:
	rm -f Test $(tests_ddMd_potentials_OBJS)
	rm -f $(tests_ddMd_potentials_OBJS:.o=.d)

clean-deps:
	rm -f $(tests_ddMd_potentials_OBJS:.o=.d)

Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(tests_ddMd_potentials_OBJS:.o=.d)

