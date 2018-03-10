BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_boundary_OBJS)

clean:
	rm -f $(simp_boundary_OBJS) $(simp_boundary_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_boundary_OBJS:.o=.d)

-include $(simp_boundary_OBJS:.o=.d)

