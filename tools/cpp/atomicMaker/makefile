
SRC_DIR_REL=../../../src
SRC_DIR=../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: AtomicMaker.o AtomicMaker

clean:	
	-rm -f AtomicMaker.o
	-rm -f AtomicMaker
	-rm -f AtomicMaker.d

-include AtomicMaker.d

AtomicMaker: AtomicMaker.cpp $(inter_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o AtomicMaker AtomicMaker.o \
               $(inter_OBJS) $(util_OBJS)
