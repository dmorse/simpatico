#-----------------------------------------------------------------------
# The copy of this namespace-level makefile in the src/ directory is 
# copied to the obj/serial and obj/parallel directories by the setup
# script to create the copies in those directories. Only the copy in
# the src/ directory is stored in the repository.
#-----------------------------------------------------------------------
# Makefile includes 

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/mdPp/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdPp/sources.mk
#-----------------------------------------------------------------------
# Main targets

#all: $(mdPp_OBJS) $(mdPp_LIB) $(MDPP).o $(MDPP_EXE)
all: $(mdPp_OBJS) $(mdPp_LIB)

clean:
	rm -f $(mdPp_OBJS) $(mdPp_OBJS:.o=.d) $(mdPp_LIB)
	rm -f $(MDPP).o $(MDPP).d

veryclean:
	make clean
	rm -f lib*.a
#ifeq ($(OBJ_DIR),$(SRC_DIR))
#	rm -f configIos/ConfigIoFactory.cpp
#	rm -f analyzers/AnalyzerFactory.cpp
#endif

#$(MDPP_EXE): $(MDPP).o $(LIBS)
#	$(CXX) $(LDFLAGS) -o $(MDPP_EXE) $(MDPP).o $(LIBS)

#mdPp:
#	make $(MDPP_EXE)

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(mdPp_OBJS:.o=.d)
#-include $(MDPP).d

