SRC_DIR=../../../src

include $(SRC_DIR)/makefile.all
include $(SRC_DIR)/mcMd/makefile.inc
include $(SRC_DIR)/util/makefile.inc

clean:	
	-rm ChainMaker

# ==========================================================================

CHAINMAKER=ChainMaker

$(CHAINMAKER): $(CHAINMAKER).cpp $(mcMd_OBJS) $(util_OBJS)
	$(CXX) $(LDFLAGS) $(DEFINES) $(INCLUDES) -o $(CHAINMAKER) $(CHAINMAKER).cpp \
                          $(mcMd_OBJS) $(util_OBJS)

all: $(CHAINMAKER)

