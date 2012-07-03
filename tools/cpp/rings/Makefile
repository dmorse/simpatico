
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

all: RingMaker.o RingMaker

clean:	
	-rm -f RingMaker.o
	-rm -f RingMaker.d
	-rm -f RingMaker

-include RingMaker.d

RingMaker: RingMaker.cpp $(inter_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o RingMaker RingMaker.o \
               $(inter_OBJS) $(util_OBJS)
