# Relative path to the simpatico/src directory
SRC_DIR_REL =../..

# Makefile includes 
include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
#-----------------------------------------------------------------------

MD_PROCESS_OBJ=$(BLD_DIR)/ddMd/sp/mdProcess.o
MD_PROCESS_EXE=$(BIN_DIR)/mdProcess

# Main targets

all: $(ddMd_sp_OBJS)

clean:
	rm -f $(ddMd_sp_OBJS) $(ddMd_sp_OBJS:.o=.d) 

clean-deps:
	rm -f $(ddMd_sp_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files
 
-include $(util_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)

