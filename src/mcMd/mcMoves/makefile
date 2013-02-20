SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
#include $(SRC_DIR_REL)/mcMd/mcMoves/sources.mk

include sources.mk

all: $(mcMd_mcMoves_OBJS)

clean:
	rm -f $(mcMd_mcMoves_OBJS) $(mcMd_mcMoves_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_mcMoves_OBJS:.o=.d)

-include $(mcMd_mcMoves_OBJS:.o=.d)

