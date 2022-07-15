BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_processor_OBJS)

clean:
	rm -f $(mdPp_processor_OBJS) $(mdPp_processor_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_processor_OBJS:.o=.d)

-include $(mdPp_processor_OBJS:.o=.d)

