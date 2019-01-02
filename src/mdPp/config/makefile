BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_config_OBJS)

clean:
	rm -f $(mdPp_config_OBJS) $(mdPp_config_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_config_OBJS:.o=.d)

-include $(mdPp_config_OBJS:.o=.d)

