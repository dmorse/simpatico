BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_commands_OBJS)

clean:
	rm -f $(mcMd_commands_OBJS) $(mcMd_commands_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_commands_OBJS:.o=.d)

-include $(mcMd_commands_OBJS:.o=.d)

