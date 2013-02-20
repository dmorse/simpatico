SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/compiler.mk
include $(SRC_DIR_REL)/util/defines.mk
include $(SRC_DIR_REL)/inter/defines.mk
include $(SRC_DIR_REL)/mcMd/defines.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/mcMd/diagnostics/util/sources.mk

all: $(mcMd_diagnostics_util_OBJS)

clean:
	rm -f $(mcMd_diagnostics_util_OBJS) $(mcMd_diagnostics_util_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_diagnostics_util_OBJS:.o=.d)

-include $(mcMd_diagnostics_util_OBJS:.o=.d)

