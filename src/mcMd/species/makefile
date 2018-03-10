BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_species_OBJS)

clean:
	rm -f $(mcMd_species_OBJS) $(mcMd_species_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_species_OBJS:.o=.d)

-include $(mcMd_species_OBJS:.o=.d)

