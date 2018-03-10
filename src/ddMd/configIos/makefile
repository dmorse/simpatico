BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_configIos_OBJS)

clean:
	rm -f $(ddMd_configIos_OBJS) $(ddMd_configIos_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_configIos_OBJS:.o=.d)

-include $(ddMd_configIos_OBJS:.o=.d)

