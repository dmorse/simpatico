BLD_DIR_REL=../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/modifiers/sources.mk

all: $(ddMd_tests_modifiers_OBJS)

clean:
	rm -f $(ddMd_tests_modifiers_OBJS) 
	rm -f $(ddMd_tests_modifiers_OBJS:.o=.d)
	rm -f $(ddMd_tests_modifiers_OBJS:.o=)

clean-deps:
	rm -f $(ddMd_tests_modifiers_OBJS:.o=.d)

-include $(ddMd_tests_modifiers_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

