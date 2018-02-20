BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_chemistry_OBJS)

clean:
	rm -f $(ddMd_chemistry_OBJS) $(ddMd_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_chemistry_OBJS:.o=.d)

-include $(ddMd_chemistry_OBJS:.o=.d)

