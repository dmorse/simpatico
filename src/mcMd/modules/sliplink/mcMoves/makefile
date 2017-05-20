# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Directory containing Simpatico library source files
BLD_DIR_REL=../../../..

#-----------------------------------------------------------------------
# Include master makefiles

include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
#include $(BLD_DIR)/util/sources.mk
#include $(BLD_DIR)/simp/sources.mk
#include $(BLD_DIR)/mcMd/sources.mk
include $(SRC_DIR)/mcMd/modules/sliplink/mcMoves/sources.mk

#-------------------------------------------------------------------
# Major targets

all: $(mcMd_modules_sliplink_mcMoves_OBJS)

clean:
	-rm -f $(mcMd_modules_sliplink_mcMoves_OBJS)
	-rm -f $(mcMd_modules_sliplink_mcMoves_OBJS:.o=.d)

-include $(mcMd_modules_sliplink_mcMoves_OBJS:.o=.d)

#-------------------------------------------------------------------
