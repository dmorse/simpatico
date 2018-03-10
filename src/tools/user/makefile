BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/tools/patterns.mk

all: $(tools_user_OBJS)

clean:
	rm -f $(tools_user_OBJS) $(tools_user_OBJS:.o=.d)

clean-deps:
	rm -f $(tools_user_OBJS:.o=.d)

-include $(tools_user_OBJS:.o=.d)

