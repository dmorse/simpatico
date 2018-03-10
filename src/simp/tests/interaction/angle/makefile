BLD_DIR_REL =../../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/tests/interaction/angle/sources.mk

all: $(simp_tests_interaction_angle_EXES)

clean:
	rm -f $(simp_tests_interaction_angle_EXES) 
	rm -f $(simp_tests_interaction_angle_OBJS)
	rm -f $(simp_tests_interaction_angle_OBJS:.o=.d)

-include $(simp_tests_interaction_angle_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

