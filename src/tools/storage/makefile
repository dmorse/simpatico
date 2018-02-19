BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_storage_OBJS)

clean:
	rm -f $(tools_storage_OBJS) $(tools_storage_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_storage_OBJS:.o=.d)

-include $(tools_storage_OBJS:.o=.d)

