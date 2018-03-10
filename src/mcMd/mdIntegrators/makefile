BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_mdIntegrators_OBJS)

clean:
	rm -f $(mcMd_mdIntegrators_OBJS) $(mcMd_mdIntegrators_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_mdIntegrators_OBJS:.o=.d)

-include $(mcMd_mdIntegrators_OBJS:.o=.d)

