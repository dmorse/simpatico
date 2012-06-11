SRC_DIR_REL =../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

#-----------------------------------------------------------------------
# Variable definitions

#-----------------------------------------------------------------------


TESTS_DIR=$(ROOT_DIR)/tests

all:  Test

clean:
	rm -f Test Test.o Test.d
	rm -f */*.o */*/*.o 
	rm -f */*.d */*/*.d

Test: $(mcMd_LIB) $(inter_LIB) $(util_LIB) Test.o
	$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.o \
	       $(mcMd_LIB) $(inter_LIB) $(util_LIB)

-include Test.d

