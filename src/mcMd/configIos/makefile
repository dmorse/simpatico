BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_configIos_OBJS)

clean:
	rm -f $(mcMd_configIos_OBJS) $(mcMd_configIos_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_configIos_OBJS:.o=.d)

-include $(mcMd_configIos_OBJS:.o=.d)

