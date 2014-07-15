BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/pair/sources.mk

all: $(inter_tests_pair_OBJS) 

clean:
	rm -f $(inter_tests_pair_OBJS)
	rm -f $(inter_tests_pair_OBJS:.o=.d)
	rm -f $(inter_tests_pair_OBJS:.o=)
	$(MAKE) clean-outputs

clean-outputs:
	rm -f out/*

-include $(inter_tests_pair_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

