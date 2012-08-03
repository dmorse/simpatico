SRC_DIR_REL =../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(TESTS_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk

all:  Test

clean:
	rm -f Test Test.o Test.d
	cd accumulators; make clean
	cd archives; make clean
	cd boundary; make clean
	cd containers; make clean
	cd crystal; make clean
	cd format; make clean
	cd param; make clean
	cd random; make clean
	cd space; make clean
	cd util; make clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

Test: $(util_LIB) $(TESTS_DIR)/util/Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test $(TESTS_DIR)/util/Test.o \
	       $(util_LIB)

-include Test.d

