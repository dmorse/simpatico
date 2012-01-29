SRC_DIR =../..

include $(SRC_DIR)/makefile.all
include $(SRC_DIR)/util/makefile.inc
include $(SRC_DIR)/mcMd/makefile.inc
include makefile.inc

.PHONY: all clean

all: $(draft_mcMd_OBJS)

clean:
	rm -f $(draft_mcMd_OBJS)

-include $(draft_mcMd_OBJS:.o=.d)

