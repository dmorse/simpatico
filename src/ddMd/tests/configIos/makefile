OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/ddMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/configIos/sources.mk

all: $(ddMd_tests_configIos_OBJS) $(ddMd_tests_configIos_OBJS:.o=)

clean:
	rm -f $(ddMd_tests_configIos_OBJS) $(ddMd_tests_configIos_OBJS:.o=.d)
	rm -f $(ddMd_tests_configIos_OBJS:.o=)
	rm -f binary out

-include $(ddMd_tests_configIos_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

