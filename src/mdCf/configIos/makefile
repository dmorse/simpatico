SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/mdCf/config.mk
include $(SRC_DIR)/mdCf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdCf/sources.mk

all: $(mdCf_configIos_OBJS)

clean:
	rm -f $(mdCf_configIos_OBJS) $(mdCf_configIos_OBJS:.o=.d)

clean-deps:
	rm -f $(mdCf_configIos_OBJS:.o=.d)

-include $(mdCf_configIos_OBJS:.o=.d)

