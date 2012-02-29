
SRC_DIR_REL=../../../src
SRC_DIR=../../../src

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/sources.mk

CHAINMAKER=ChainMaker

all: $(CHAINMAKER).o $(CHAINMAKER)

clean:	
	-rm ChainMaker.o
	-rm ChainMaker.d
	-rm ChainMaker

-include $(CHAINMAKER).d

$(CHAINMAKER): $(CHAINMAKER).cpp $(mcMd_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $(CHAINMAKER) $(CHAINMAKER).o \
               $(mcMd_OBJS) $(util_OBJS)
