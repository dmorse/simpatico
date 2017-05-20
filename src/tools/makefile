#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes 

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/tools/sources.mk

#-----------------------------------------------------------------------
# Executable program targets

# MD postprocessor --------------------------------
mdPp=$(BLD_DIR)/tools/mdPp

$(mdPp_BIN): $(mdPp).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(mdPp_BIN) $(mdPp).o $(LIBS)

mdPp:
	$(MAKE) $(mdPp_BIN)

-include $(mdPp).d

# ChainMaker --------------------------------------

ChainMaker=$(BLD_DIR)/tools/generators/ChainMaker
ChainMaker_BIN=$(BIN_DIR)/chainMaker

$(ChainMaker_BIN): $(ChainMaker).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(ChainMaker_BIN) $(ChainMaker).o $(LIBS)

-include $(ChainMaker).d

#-----------------------------------------------------------------------
# Main targets

all: $(tools_LIB) $(mdPp_BIN) $(ChainMaker_BIN)

clean:
	rm -f $(tools_OBJS) $(tools_OBJS:.o=.d) $(tools_LIB)
	rm -f $(mdPp).o $(mdPp).d
	rm -f $(ChainMaker).o $(ChainMaker).d
	cd tests; $(MAKE) clean
	rm -f *.o */*.o 
	rm -f *.d */*.d 

veryclean:
	$(MAKE) clean
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f config/ConfigReaderFactory.cpp
	rm -f config/ConfigWriterFactory.cpp
	rm -f trajectory/TrajectoryReaderFactory.cpp
	rm -f processor/ProcessorAnalyzerFactory.cpp
endif

#-----------------------------------------------------------------------
# Include class dependency files
 
-include $(util_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)

