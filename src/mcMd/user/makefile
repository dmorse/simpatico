BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/mcMd/patterns.mk

all: $(mcMd_user_OBJS)

clean:
	rm -f $(mcMd_user_OBJS) $(mcMd_user_OBJS:.o=.d)

clean-deps:
	rm -f $(mcMd_user_OBJS:.o=.d)

-include $(mcMd_user_OBJS:.o=.d)

