BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/simp/patterns.mk
include $(SRC_DIR)/simp/user/sources.mk

all: $(simp_user_OBJS)

clean:
	rm -f $(simp_user_OBJS) $(simp_user_OBJS:.o=.d)

clean-deps:
	rm -f $(simp_user_OBJS:.o=.d)

-include $(simp_user_OBJS:.o=.d)

