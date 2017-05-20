# ------------------------------------------------------------------
# Users may need to change the following Makefile macros:
#  - SRC_DIR_REL, if this directory is moved or copied to another location.
#  - BIN_DIR, to change the directory for executables

# Path to directory containing Simpatico library source files
BLD_DIR_REL=../../..

# Include master makefiles
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
#include $(SRC_DIR)/util/sources.mk
#include $(SRC_DIR)/simp/sources.mk
#include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/mcMd/modules/sliplink/sources.mk

#-------------------------------------------------------------------
# Major targets

.phony: all clean

all: $(mcMd_modules_sliplink_OBJS)

clean:
	-rm -f $(mcMd_modules_sliplink_OBJS) $(mcMd_modules_sliplink_OBJS:.o=.d)

-include $(modules_sliplink_OBJS:.o=.d)
