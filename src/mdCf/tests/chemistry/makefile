OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(OBJ_DIR)/mdCf/config.mk
include $(SRC_DIR)/mdCf/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mdCf/sources.mk
include $(SRC_DIR)/mdCf/tests/chemistry/sources.mk

all: $(ddMd_tests_sp_chemistry_OBJS)

clean:
	rm -f $(ddMd_tests_sp_chemistry_OBJS) 
	rm -f $(ddMd_tests_sp_chemistry_OBJS:.o=.d)
	rm -f $(ddMd_tests_sp_chemistry_OBJS:.o=)

-include $(ddMd_tests_sp_chemistry_OBJS:.o=.d)
-include $(ddMd_sp_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

