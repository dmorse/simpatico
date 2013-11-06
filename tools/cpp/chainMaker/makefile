
SRC_DIR_REL=../../../src
SRC_DIR=../../../src

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/inter/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/inter/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

all: ChainMaker.o ChainMaker

clean:	
	-rm -f ChainMaker.o
	-rm -f ChainMaker.d
	-rm -f ChainMaker

-include ChainMaker.d

ChainMaker: ChainMaker.cpp $(inter_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o ChainMaker ChainMaker.o \
               $(inter_OBJS) $(util_OBJS)
