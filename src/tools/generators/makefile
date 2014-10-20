SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/tools/sources.mk

PROGS=$(BIN_DIR)/ChainMaker\
      $(BIN_DIR)/AtomicMaker\
      $(BIN_DIR)/RingMaker

all: $(tools_generators_OBJS) $(PROGS)

clean:
	rm -f $(tools_generators_OBJS) $(tools_generators_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_generators_OBJS:.o=.d)

ChainMaker=$(BLD_DIR)/tools/generators/chainMaker/ChainMaker
$(BIN_DIR)/ChainMaker: $(ChainMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(BIN_DIR)/ChainMaker $(ChainMaker).o $(LIBS)

AtomicMaker=$(BLD_DIR)/tools/generators/atomicMaker/AtomicMaker
$(BIN_DIR)/AtomicMaker: $(AtomicMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(BIN_DIR)/AtomicMaker $(AtomicMaker).o $(LIBS)

RingMaker=$(BLD_DIR)/tools/generators/rings/RingMaker
$(BIN_DIR)/RingMaker: $(RingMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(BIN_DIR)/RingMaker $(RingMaker).o $(LIBS)

-include $(tools_generators_OBJS:.o=.d)

