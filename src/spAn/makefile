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
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/spAn/sources.mk

#-----------------------------------------------------------------------
# Variable definition

# MD postprocessor program location
MDPP=$(BLD_DIR)/spAn/mdPp

#-----------------------------------------------------------------------
# Main targets

all: $(spAn_OBJS) $(spAn_LIB) $(MDPP).o $(MDPP_EXE)

$(MDPP_EXE): $(MDPP).o $(LIBS)
	$(CXX) $(LDFLAGS) -o $(MDPP_EXE) $(MDPP).o $(LIBS)

mdPp:
	$(MAKE) $(MDPP_EXE)

clean:
	rm -f $(spAn_OBJS) $(spAn_OBJS:.o=.d) $(spAn_LIB)
	rm -f $(MDPP).o $(MDPP).d
	cd tests; $(MAKE) clean
	rm -f *.o */*.o 
	rm -f *.d */*.d 

veryclean:
	$(MAKE) clean
	rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/AnalyzerFactory.cpp
endif

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(spAn_OBJS:.o=.d)
-include $(MDPP).d

