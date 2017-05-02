SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include sources.mk

all: $(mcMd_species_OBJS)

clean:
	rm -f $(mcMd_species_OBJS) $(mcMd_species_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_species_OBJS:.o=.d)

-include $(mcMd_species_OBJS:.o=.d)

