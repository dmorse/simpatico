# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Directory containing Simpatico library source files
BLD_DIR_REL=../../..

# Parent directory for modules
MODULES_DIR=$(BLD_DIR)/modules

#-----------------------------------------------------------------------
# Include master makefiles

include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(MODULES_DIR)/sliplink/analyzers/sources.mk

#-------------------------------------------------------------------
# Major targets

all: $(modules_sliplink_analyzers_OBJS)

clean:
	-rm -f $(modules_sliplink_analyzers_OBJS)
	-rm -f $(modules_sliplink_analyzers_OBJS:.o=.d)

-include $(modules_sliplink_analyzers_OBJS:.o=.d)

#-------------------------------------------------------------------
