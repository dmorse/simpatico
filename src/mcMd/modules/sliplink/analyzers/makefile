# Directory containing Simpatico library source files
BLD_DIR_REL=../../../..

#-----------------------------------------------------------------------
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

all: $(mcMd_modules_sliplink_analyzers_OBJS)

clean:
	-rm -f $(mcMd_modules_sliplink_analyzers_OBJS)
	-rm -f $(mcMd_modules_sliplink_analyzers_OBJS:.o=.d)

-include $(mcMd_modules_sliplink_analyzers_OBJS:.o=.d)

#-------------------------------------------------------------------
