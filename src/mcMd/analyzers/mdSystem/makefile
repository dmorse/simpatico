SRC_DIR_REL =../../..

include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR_REL)/util/config.mk
include $(SRC_DIR_REL)/simp/config.mk
include $(SRC_DIR_REL)/mcMd/config.mk
include $(SRC_DIR_REL)/mcMd/patterns.mk
include $(SRC_DIR_REL)/util/sources.mk
include $(SRC_DIR_REL)/mcMd/analyzers/mdSystem/sources.mk

all: $(mcMd_analyzers_mdSystem_OBJS)

clean:
	rm -f $(mcMd_analyzers_mdSystem_OBJS) $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

-include $(mcMd_analyzers_mdSystem_OBJS:.o=.d)

