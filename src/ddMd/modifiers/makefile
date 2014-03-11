OBJ_DIR_REL =../..

include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/patterns.mk
include $(OBJ_DIR)/inter/patterns.mk
include $(OBJ_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/ddMd/modifiers/sources.mk

all: $(ddMd_modifiers_OBJS)

clean:
	rm -f $(ddMd_modifiers_OBJS) $(ddMd_modifiers_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_modifiers_OBJS:.o=.d)

-include $(ddMd_modifiers_OBJS:.o=.d)

