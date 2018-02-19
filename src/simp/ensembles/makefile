BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_ensembles_OBJS)

clean:
	rm -f $(simp_ensembles_OBJS) $(simp_ensembles_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_ensembles_OBJS:.o=.d)

-include $(simp_ensembles_OBJS:.o=.d)

