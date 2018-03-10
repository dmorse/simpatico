BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_processor_OBJS)

clean:
	rm -f $(tools_processor_OBJS) $(tools_processor_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_processor_OBJS:.o=.d)

-include $(tools_processor_OBJS:.o=.d)

