BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_generators_OBJS)

clean:
	rm -f $(mcMd_generators_OBJS) $(mcMd_generators_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_generators_OBJS:.o=.d)

-include $(mcMd_generators_OBJS:.o=.d)

