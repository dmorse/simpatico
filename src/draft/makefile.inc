# Note: The SRC macro must be defined by any Makefile that includes this file.

include $(SRC)/mcMd/makefile.inc

# Subdirectory macros
DRAFT_MCMD_DIR=$(SRC)/draft/mcMd

# Include subdirectory makefile.local files
include $(DRAFT_MCMD_DIR)/makefile.local

HDR_MCMD_ALL+=$(DRAFT_MCMD_HDR)
SRC_MCMD_ALL+=$(DRAFT_MCMD_SRC)
OBJ_MCMD_ALL+=$(DRAFT_MCMD_OBJ)
