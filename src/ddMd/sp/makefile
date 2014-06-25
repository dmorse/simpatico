# Relative path to the simpatico/src directory
SRC_DIR_REL =../..

# Makefile includes 
include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
#-----------------------------------------------------------------------

MD_PROCESS_OBJ=$(OBJ_DIR)/ddMd/sp/mdProcess.o
MD_PROCESS_EXE=$(BIN_DIR)/mdProcess

# Main targets

all: $(ddMd_sp_OBJS)

clean:
	rm -f $(ddMd_sp_OBJS) $(ddMd_sp_OBJS:.o=.d) 

veryclean:
	make clean
ifeq ($(OBJ_DIR),$(SRC_DIR))
	rm -f configIos/ConfigIoFactory.cpp
	rm -f analyzers/AnalyzerFactory.cpp
endif

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_sp_OBJS:.o=.d)

