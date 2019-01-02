BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

# -- Targets for programs

ChainMaker=$(BLD_DIR)/mdPp/generators/ChainMaker
ChainMaker_BIN=$(BIN_DIR)/chainMaker

$(ChainMaker_BIN): $(ChainMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(ChainMaker_BIN) $(ChainMaker).o $(LIBS)

all: $(ChainMaker_BIN) 

clean:
	rm -f $(ChainMaker).o $(ChainMaker).d

clean-deps:
	rm -f $(ChainMaker).d

-include $(mdPp_OBJS:.o=.d)
-include $(ChainMaker).d

