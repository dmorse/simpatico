#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the obj/serial and obj/parallel directories by the setup script.
# Only the copy in the src/ directory is stored in the repository. 
#-----------------------------------------------------------------------
# Relative path to the simpatico/src directory
SRC_DIR_REL =..

# Makefile includes 
include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/mdCf/config.mk
include $(SRC_DIR)/mdCf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdCf/sources.mk
#-----------------------------------------------------------------------
# Main targets

all: $(ddMd_sp_OBJS) $(ddMd_sp_LIB)
#all: $(ddMd_sp_OBJS) $(ddMd_sp_LIB) $(MDPP).o $(MDPP_EXE)

clean:
	rm -f $(ddMd_sp_OBJS) $(ddMd_sp_OBJS:.o=.d) $(ddMd_sp_LIB)
#	rm -f $(MDPP).o $(MDPP).d

veryclean:
	make clean
	rm -f lib*.a
ifeq ($(OBJ_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/AnalyzerFactory.cpp
endif

#$(MDPP_EXE): $(MDPP).o $(LIBS)
#	$(CXX) $(LDFLAGS) -o $(MDPP_EXE) $(MDPP).o $(LIBS)

#mdPp:
#	make $(MDPP_EXE)

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_sp_OBJS:.o=.d)
#-include $(MDPP).d

