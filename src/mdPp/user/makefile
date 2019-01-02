BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mdPp/patterns.mk

all: $(mdPp_user_OBJS)

clean:
	rm -f $(mdPp_user_OBJS) $(mdPp_user_OBJS:.o=.d)

clean-deps:
	rm -f $(mdPp_user_OBJS:.o=.d)

-include $(mdPp_user_OBJS:.o=.d)

