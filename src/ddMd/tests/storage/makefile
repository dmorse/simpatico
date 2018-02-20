BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/tests/storage/sources.mk

all: $(ddMd_tests_storage_EXES)

clean:
	rm -f $(ddMd_tests_storage_EXES)
	rm -f $(ddMd_tests_storage_OBJS) 
	rm -f $(ddMd_tests_storage_OBJS:.o=.d)

-include $(ddMd_tests_storage_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

