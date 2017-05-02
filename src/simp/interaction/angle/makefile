SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/simp/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/interaction/angle/sources.mk

all: $(simp_interaction_angle_OBJS)

clean:
	rm -f $(simp_interaction_angle_OBJS) $(simp_interaction_angle_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_interaction_angle_OBJS:.o=.d)

-include $(simp_interaction_angle_OBJS:.o=.d)

