SRC_DIR_REL =../../src
include $(SRC_DIR_REL)/compiler.mk

include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/mcMd/defines.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

TESTS_DIR=$(ROOT_DIR)/tests
include $(TESTS_DIR)/mcMd/patterns.mk

#-----------------------------------------------------------------------
# Variable definitions

#-----------------------------------------------------------------------




all:  Test

clean:
	-rm -f Test Test.o Test.d
	-rm -f */*.o */*/*.o 
	-rm -f */*.d */*/*.d
	-rm -f *.cfg
	-rm -f log

Test: $(mcMd_LIB) $(inter_LIB) $(util_LIB) Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include Test.d

