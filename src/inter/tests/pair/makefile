OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/pair/sources.mk

TEST=inter/tests/pair/Test

all: $(inter_tests_pair_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(inter_tests_pair_OBJS) $(inter_tests_pair_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 
	make clean-outputs

clean-outputs:
	rm -f out/*

-include $(inter_tests_pair_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

