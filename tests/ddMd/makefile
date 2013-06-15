SRC_DIR_REL =../../src
include $(SRC_DIR_REL)/compiler.mk

include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/ddMd/defines.mk
include $(TESTS_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

#-----------------------------------------------------------------------


all:  Test

clean:
	cd chemistry; make clean
	cd storage; make clean
	cd communicate; make clean
	cd neighbor; make clean
	cd potentials; make clean
	cd configIos; make clean
	cd simulation; make clean
	-rm -f Test.o Test.d Test
	-rm -f log

Test: $(ddMd_LIB) $(inter_LIB) $(util_LIB) $(TESTS_DIR)/ddMd/Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test $(TESTS_DIR)/ddMd/Test.o \
	       $(ddMd_LIB) $(inter_LIB) $(util_LIB)

-include Test.d

