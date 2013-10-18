BLD_DIR_REL =../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/ddMd/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/sources.mk

all: $(ddMd_tests_OBJS) Test

clean:
	rm -f $(ddMd_tests_OBJS) Test $(ddMd_tests_OBJS:.o=.d)
	cd chemistry; make clean
	cd storage; make clean
	cd communicate; make clean
	cd neighbor; make clean
	cd potentials; make clean
	cd configIos; make clean
	cd simulation; make clean
	-rm -f Test
	-rm -f log


clean-deps:
	rm -f $(ddMd_tests_OBJS:.o=.d)
	*/*.d

Test: Test.o $(ddMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include $(ddMd_tests_OBJS:.o=.d)


