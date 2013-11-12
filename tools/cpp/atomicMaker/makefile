
SRC_DIR_REL=../../../src
SRC_DIR=../../../src

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/util/config.mk
include $(SRC_DIR)/inter/config.mk
include $(SRC_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mcMd/sources.mk

all: AtomicMaker.o AtomicMaker

clean:	
	-rm -f AtomicMaker
	-rm -f AtomicMaker.o
	-rm -f AtomicMaker.d

-include AtomicMaker.d

AtomicMaker: AtomicMaker.cpp $(inter_LIB) $(util_LIB)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o AtomicMaker AtomicMaker.o \
               $(inter_LIB) $(util_LIB)
