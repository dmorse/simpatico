BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_analyzers_OBJS)

clean:
	rm -f $(tools_analyzers_OBJS) $(tools_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_analyzers_OBJS:.o=.d)

-include $(tools_analyzers_OBJS:.o=.d)

