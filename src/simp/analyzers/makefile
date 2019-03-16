BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_analyzers_OBJS)

clean:
	rm -f $(simp_analyzers_OBJS) $(simp_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_analyzers_OBJS:.o=.d)

-include $(simp_analyzers_OBJS:.o=.d)

