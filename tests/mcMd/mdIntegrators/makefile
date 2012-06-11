SRC_DIR_REL=../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all:    Test

clean:	 
	rm -f Test 

Test: Test.cpp MdIntegratorTest.h
	$(CXX) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.cpp  $(util_OBJS) $(mcMd_OBJS)

