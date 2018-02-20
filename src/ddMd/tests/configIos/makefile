BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/configIos/sources.mk

all: $(ddMd_tests_configIos_EXES) 

clean:
	rm -f $(ddMd_tests_configIos_EXES)
	rm -f $(ddMd_tests_configIos_OBJS) 
	rm -f $(ddMd_tests_configIos_OBJS:.o=.d)
	rm -f tmp/*

clean-outputs:
	@rm -f tmp/*

-include $(ddMd_tests_configIos_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

