#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the bld/serial and bld/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes 

BLD_DIR_REL =..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

#-----------------------------------------------------------------------
# Main targets

all: $(tools_LIB) $(mdPp_EXE) $(chainMaker_EXE)

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
# Executable program targets

#--- Postprocessor --------------------------------
# Note: mdPp_EXE is defined in tools/config.mk

mdPp=$(BLD_DIR)/tools/mdPp

$(mdPp_EXE): $(mdPp).o $(LIBS)
	$(CXX) -o $(mdPp_EXE) $(mdPp).o $(LIBS) $(LDFLAGS)

mdPp:
	$(MAKE) $(mdPp_EXE)


#--- ChainMaker -------------------------------------
# Note: chainMaker_EXE is defined in tools/config.mk

ChainMaker=$(BLD_DIR)/tools/generators/ChainMaker

$(chainMaker_EXE): $(ChainMaker).o $(LIBS)
	$(CXX) -o $(chainMaker_EXE) $(ChainMaker).o $(LIBS) $(LDFLAGS)

#-----------------------------------------------------------------------
# Include *.d dependency files
 
-include $(util_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(mdPp).d
-include $(ChainMaker).d

