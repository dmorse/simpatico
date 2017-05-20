SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/simp/patterns.mk
include $(SRC_DIR_REL)/simp/interaction/pair/sources.mk

all: $(simp_interaction_pair_OBJS)

clean:
	rm -f $(simp_interaction_pair_OBJS) $(simp_interaction_pair_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_interaction_pair_OBJS:.o=.d)

-include $(simp_interaction_pair_OBJS:.o=.d)

