SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/mcMd/config.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/analyzers/sources.mk

all: $(mcMd_analyzers_OBJS)

clean:
	rm -f $(mcMd_analyzers_OBJS) $(mcMd_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_analyzers_OBJS:.o=.d)

-include $(mcMd_analyzers_OBJS:.o=.d)

