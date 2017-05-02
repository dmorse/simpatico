
SRC_DIR_REL=../../../src
SRC_DIR=../../../src

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/simp/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: RingMaker.o RingMaker

clean:	
	-rm -f RingMaker.o
	-rm -f RingMaker.d
	-rm -f RingMaker

-include RingMaker.d

RingMaker: RingMaker.cpp $(simp_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o RingMaker RingMaker.o \
               $(simp_OBJS) $(util_OBJS)
