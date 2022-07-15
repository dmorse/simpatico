BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_crystal_OBJS)

clean:
	rm -f $(simp_crystal_OBJS) $(simp_crystal_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_crystal_OBJS:.o=.d)

-include $(simp_crystal_OBJS:.o=.d)

