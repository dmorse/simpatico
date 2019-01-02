BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_analyzers_OBJS)

clean:
	rm -f $(mdPp_analyzers_OBJS) $(mdPp_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_analyzers_OBJS:.o=.d)

-include $(mdPp_analyzers_OBJS:.o=.d)

