SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/mdPp/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/mdPp/sources.mk

all: $(mdPp_chemistry_OBJS)

clean:
	rm -f $(mdPp_chemistry_OBJS) $(mdPp_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_chemistry_OBJS:.o=.d)

-include $(mdPp_chemistry_OBJS:.o=.d)

