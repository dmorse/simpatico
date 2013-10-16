SRC_DIR_REL =../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR)/util/defines.mk
include $(SRC_DIR)/inter/defines.mk
include $(SRC_DIR)/mcMd/defines.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/diagnostics/sources.mk

all: $(mcMd_diagnostics_OBJS)

clean:
	rm -f $(mcMd_diagnostics_OBJS) $(mcMd_diagnostics_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_diagnostics_OBJS:.o=.d)

-include $(mcMd_diagnostics_OBJS:.o=.d)

