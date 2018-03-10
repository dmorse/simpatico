BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_misc_OBJS)

clean:
	rm -f $(ddMd_misc_OBJS) $(ddMd_misc_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_misc_OBJS:.o=.d)

-include $(ddMd_misc_OBJS:.o=.d)

