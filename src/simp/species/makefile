BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk

all: $(simp_species_OBJS)

clean:
	rm -f $(simp_species_OBJS) $(simp_species_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_species_OBJS:.o=.d)

-include $(simp_species_OBJS:.o=.d)

