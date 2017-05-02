SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/interaction/external/sources.mk

all: $(simp_interaction_external_OBJS)

clean:
	rm -f $(simp_interaction_external_OBJS) $(simp_interaction_external_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_interaction_external_OBJS:.o=.d)

-include $(simp_interaction_external_OBJS:.o=.d)

