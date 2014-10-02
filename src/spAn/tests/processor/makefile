BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/spAn/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/spAn/tests/processor/sources.mk

all: $(spAn_tests_processor_OBJS)

clean:
	rm -f $(spAn_tests_processor_OBJS) 
	rm -f $(spAn_tests_processor_OBJS:.o=.d)
	rm -f $(spAn_tests_processor_OBJS:.o=)
	rm -f out/*

-include $(spAn_tests_processor_OBJS:.o=.d)
-include $(spAn_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

