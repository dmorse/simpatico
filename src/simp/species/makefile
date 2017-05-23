SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/simp/patterns.mk
include sources.mk

all: $(simp_species_OBJS)

clean:
	rm -f $(simp_species_OBJS) $(simp_species_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_species_OBJS:.o=.d)

-include $(simp_species_OBJS:.o=.d)

