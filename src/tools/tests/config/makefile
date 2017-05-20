BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/tools/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/tools/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/tools/tests/config/sources.mk

all: $(tools_tests_config_OBJS)

clean:
	rm -f $(tools_tests_config_OBJS) 
	rm -f $(tools_tests_config_OBJS:.o=.d)
	rm -f $(tools_tests_config_OBJS:.o=)
	rm -f out/*

-include $(tools_tests_config_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

