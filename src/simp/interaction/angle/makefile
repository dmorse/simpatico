BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_interaction_angle_OBJS)

clean:
	rm -f $(simp_interaction_angle_OBJS) $(simp_interaction_angle_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_interaction_angle_OBJS:.o=.d)

-include $(simp_interaction_angle_OBJS:.o=.d)

