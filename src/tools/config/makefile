BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_config_OBJS)

clean:
	rm -f $(tools_config_OBJS) $(tools_config_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_config_OBJS:.o=.d)

-include $(tools_config_OBJS:.o=.d)

