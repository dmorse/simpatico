SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdPp/sources.mk

all: $(mdPp_analyzers_OBJS)

clean:
	rm -f $(mdPp_analyzers_OBJS) $(mdPp_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_analyzers_OBJS:.o=.d)

-include $(mdPp_analyzers_OBJS:.o=.d)

