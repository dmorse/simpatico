BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_analysis_OBJS)

clean:
	rm -f $(simp_analysis_OBJS) $(simp_analysis_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_analysis_OBJS:.o=.d)

-include $(simp_analysis_OBJS:.o=.d)

