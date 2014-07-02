SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/ddMd/sources.mk

all: $(ddMd_sp_analyzers_OBJS)

clean:
	rm -f $(ddMd_sp_analyzers_OBJS) $(ddMd_sp_analyzers_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_sp_analyzers_OBJS:.o=.d)

-include $(ddMd_sp_analyzers_OBJS:.o=.d)

