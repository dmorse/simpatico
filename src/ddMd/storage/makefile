BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_storage_OBJS)

clean:
	rm -f $(ddMd_storage_OBJS) $(ddMd_storage_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_storage_OBJS:.o=.d)

-include $(ddMd_storage_OBJS:.o=.d)

