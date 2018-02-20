BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk

all: $(ddMd_integrators_OBJS)

clean:
	rm -f $(ddMd_integrators_OBJS) $(ddMd_integrators_OBJS:.o=.d)

clean-deps:
	rm -f $(ddMd_integrators_OBJS:.o=.d)

-include $(ddMd_integrators_OBJS:.o=.d)

