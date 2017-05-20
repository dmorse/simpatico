BLD_DIR_REL =../..

include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

all: $(ddMd_modifiers_OBJS)

clean:
	rm -f $(ddMd_modifiers_OBJS) $(ddMd_modifiers_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_modifiers_OBJS:.o=.d)

-include $(ddMd_modifiers_OBJS:.o=.d)

