BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_links_OBJS)

clean:
	rm -f $(mcMd_links_OBJS) $(mcMd_links_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_links_OBJS:.o=.d)

-include $(mcMd_links_OBJS:.o=.d)

