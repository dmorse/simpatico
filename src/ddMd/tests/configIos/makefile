BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/configIos/sources.mk

TEST= $(SRC_DIR)/ddMd/tests/configIos/ConfigIoTest

.Secondary: $(ddMd_tests_configIos_OBJS) 

all: $(ddMd_tests_configIos_EXES) 

run: $(TEST)
	$(MPIRUN) 6 $(TEST)

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

