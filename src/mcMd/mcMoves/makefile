SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
#include $(SRC_DIR)/mcMd/mcMoves/sources.mk

include sources.mk

all: $(mcMd_mcMoves_OBJS)

clean:
	rm -f $(mcMd_mcMoves_OBJS) $(mcMd_mcMoves_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_mcMoves_OBJS:.o=.d)

-include $(mcMd_mcMoves_OBJS:.o=.d)

