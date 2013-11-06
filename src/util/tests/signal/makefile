OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/signal/sources.mk

TEST=util/tests/signal/Test

all: $(util_tests_signal_OBJS) $(OBJ_DIR)/$(Test)

clean:
	rm -f $(util_tests_signal_OBJS) $(util_tests_signal_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 

-include $(util_tests_signal_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

