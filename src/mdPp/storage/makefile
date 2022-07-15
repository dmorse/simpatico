BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_storage_OBJS)

clean:
	rm -f $(mdPp_storage_OBJS) $(mdPp_storage_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_storage_OBJS:.o=.d)

-include $(mdPp_storage_OBJS:.o=.d)

