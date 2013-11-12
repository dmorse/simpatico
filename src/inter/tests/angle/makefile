OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/angle/sources.mk

TEST=inter/tests/angle/Test

all: $(inter_tests_angle_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(inter_tests_angle_OBJS) $(inter_tests_angle_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 

-include $(inter_tests_angle_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

