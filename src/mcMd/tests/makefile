BLD_DIR_REL =../..
include $(BLD_DIR_REL)/compiler.mk
include $(BLD_DIR)/mcMd/defines.mk
include $(BLD_DIR)/inter/defines.mk
include $(BLD_DIR)/util/defines.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/sources.mk

all: $(mcMd_tests_OBJS) Test

clean:
	rm -f $(mcMd_tests_OBJS) Test $(mcMd_tests_OBJS:.o=.d)
	cd chemistry; make clean
	cd neighbor; make clean
	cd simulation; make clean
	cd mcSimulation; make clean
	cd mdSimulation; make clean
	cd species; make clean
	-rm -f */*.o */*/*.o 
	-rm -f */*.d */*/*.d
	-rm -f *.cfg
	-rm -f log

clean-deps:
	rm -f $(mcMd_tests_OBJS:.o=.d)

Test: Test.o $(mcMd_LIB) $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include $(mcMd_tests_OBJS:.o=.d)

