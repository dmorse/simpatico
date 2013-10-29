OBJ_DIR_REL =../..

include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/ddMd/actors/sources.mk

all: $(ddMd_actors_OBJS)

clean:
	rm -f $(ddMd_actors_OBJS) $(ddMd_actors_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_actors_OBJS:.o=.d)

-include $(ddMd_actors_OBJS:.o=.d)

