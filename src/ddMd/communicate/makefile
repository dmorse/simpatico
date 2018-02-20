BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_communicate_OBJS)

clean:
	rm -f $(ddMd_communicate_OBJS) $(ddMd_communicate_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_communicate_OBJS:.o=.d)

-include $(ddMd_communicate_OBJS:.o=.d)

