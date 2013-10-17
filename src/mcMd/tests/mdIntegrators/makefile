SRC_DIR_REL=../../../src
include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/mcMd/defines.mk
include $(TESTS_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

all:    Test

clean:	 
	rm -f Test 

Test: Test.cc MdIntegratorTest.h
	$(CXX) $(TESTFLAGS) $(INCLUDES) $(DEFINES) -o Test Test.cc  \
               $(util_LIB) $(inter_LIB) $(mcMd_LIB)

