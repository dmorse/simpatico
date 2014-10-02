SRC_DIR_REL =../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/spAn/sources.mk

all: $(spAn_chemistry_OBJS)

clean:
	rm -f $(spAn_chemistry_OBJS) $(spAn_chemistry_OBJS:.o=.d)

clean-deps:
	rm -f $(spAn_chemistry_OBJS:.o=.d)

-include $(spAn_chemistry_OBJS:.o=.d)

