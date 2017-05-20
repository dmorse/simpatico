SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/tools/sources.mk

# -- Targets for programs

ChainMaker=$(BLD_DIR)/tools/generators/ChainMaker
ChainMaker_BIN=$(BIN_DIR)/chainMaker

-include $(tools_OBJS:.o=.d)
-include $(ChainMaker).d

$(ChainMaker_BIN): $(ChainMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(ChainMaker_BIN) $(ChainMaker).o $(LIBS)

all: $(ChainMaker_BIN) 

clean:
	rm -f $(ChainMaker).o $(ChainMaker).d

clean-deps:
	rm -f $(ChainMaker).d

