BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_chemistry_OBJS)

clean:
	rm -f $(mdPp_chemistry_OBJS) $(mdPp_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_chemistry_OBJS:.o=.d)

-include $(mdPp_chemistry_OBJS:.o=.d)

