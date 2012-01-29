SRC_DIR=../../..

include $(SRC_DIR)/makefile.all
include $(SRC_DIR)/util/makefile.inc
include $(SRC_DIR)/mcMd/makefile.inc
include $(SRC_DIR)/test/makefile.inc

all:    Test

clean:	 
	rm -f Test 

Test: Test.cpp MdIntegratorTest.h
	$(CXX) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.cpp  $(util_OBJS) $(mcMd_OBJS)

