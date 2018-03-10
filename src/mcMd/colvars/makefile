BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_colvars_OBJS)

clean:
	rm -f $(mcMd_colvars_OBJS) $(mcMd_colvars_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_colvars_OBJS:.o=.d)

-include $(mcMd_colvars_OBJS:.o=.d)

