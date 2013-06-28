SRC_DIR_REL =../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(TESTS_DIR)/inter/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include sources.mk

all: $(tests_inter_OBJS) Test

clean:
	rm -f $(tests_inter_OBJS) Test $(tests_inter_OBJS:.o=.d)
	cd pair; make clean;
	cd bond; make clean;
	cd angle; make clean;
	cd dihedral; make clean;
	-rm -f log

clean-deps:
	rm -f $(tests_inter_OBJS:.o=.d)

Test: Test.o $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(inter_LIB) $(util_LIB)

-include $(tests_inter_OBJS:.o=.d)

