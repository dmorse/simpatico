BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_modifiers_OBJS)

clean:
	rm -f $(ddMd_modifiers_OBJS) $(ddMd_modifiers_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_modifiers_OBJS:.o=.d)

-include $(ddMd_modifiers_OBJS:.o=.d)

