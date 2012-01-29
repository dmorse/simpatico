# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Directory containing Simpatico library source files
SRC_DIR=../../..

# Parent directory for modules
MODULES_DIR=../..

#-----------------------------------------------------------------------
# Include master makefiles

include $(SRC_DIR)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/mcMd/defines.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(MODULES_DIR)/sliplink/diagnostics/sources.mk

#-------------------------------------------------------------------
# Major targets

all: $(modules_sliplink_diagnostics_OBJS)

clean:
	-rm -f $(modules_sliplink_diagnostics_OBJS)
	-rm -f $(modules_sliplink_diagnostics_OBJS:.o=.d)

clean-deps:
	-rm -f $(modules_sliplink_diagnostics_OBJS:.o=.d)

-include $(modules_sliplink_diagnostics_OBJS:.o=.d)

#-------------------------------------------------------------------
