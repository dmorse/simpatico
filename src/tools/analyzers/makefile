SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/tools/sources.mk

all: $(tools_analyzers_OBJS)

clean:
	rm -f $(tools_analyzers_OBJS) $(tools_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_analyzers_OBJS:.o=.d)

-include $(tools_analyzers_OBJS:.o=.d)

