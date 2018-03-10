BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_chemistry_OBJS)

clean:
	rm -f $(tools_chemistry_OBJS) $(tools_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_chemistry_OBJS:.o=.d)

-include $(tools_chemistry_OBJS:.o=.d)

