BLD_DIR_REL =../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

all: $(util_tests_OBJS) Test

clean:
	rm -f $(util_tests_OBJS) Test $(util_tests_OBJS:.o=.d)
	cd accumulators/unit; make clean
	cd archives; make clean
	cd boundary; make clean
	cd containers; make clean
	cd crystal; make clean
	cd format; make clean
	cd param; make clean
	cd random; make clean
	cd space; make clean
	cd misc; make clean
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d
	-rm -f dummy
	-rm -f log

clean-deps:
	rm -f $(util_tests_OBJS:.o=.d)

Test: Test.o $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(util_LIB)

-include $(util_tests_OBJS:.o=.d)

