
SRC_DIR_REL=../../../src

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk

LIBS= $(inter_LIB) $(util_LIB)
DEFINES=$(UTIL_DEFS) $(INTER_DEFS)

# Dependencies of source files in src/inter on makefile fragments
MAKE_DEPS= -A$(BLD_DIR)/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/util/config.mk
MAKE_DEPS+= -A$(BLD_DIR)/inter/config.mk

all: ChainMaker.o ChainMaker

clean:	
	-rm -f ChainMaker.o
	-rm -f ChainMaker.d
	-rm -f ChainMaker

ChainMaker.o: ChainMaker.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c -o ChainMaker.o ChainMaker.cpp 
	
ChainMaker: ChainMaker.o $(LIBS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o ChainMaker ChainMaker.o $(LIBS)

-include ChainMaker.d
