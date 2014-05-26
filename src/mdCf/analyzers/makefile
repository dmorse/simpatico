SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(SRC_DIR)/mdCf/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdCf/sources.mk

all: $(mdCf_analyzers_OBJS)

clean:
	rm -f $(mdCf_analyzers_OBJS) $(mdCf_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(mdCf_analyzers_OBJS:.o=.d)

-include $(mdCf_analyzers_OBJS:.o=.d)

