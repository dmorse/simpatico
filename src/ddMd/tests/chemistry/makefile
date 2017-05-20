BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/ddMd/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/chemistry/sources.mk

all: $(ddMd_tests_chemistry_OBJS) $(ddMd_tests_chemistry_OBJS:.o=) 

clean:
	rm -f $(ddMd_tests_chemistry_OBJS) 
	rm -f $(ddMd_tests_chemistry_OBJS:.o=.d)
	rm -f $(ddMd_tests_chemistry_OBJS:.o=)

clean-outputs:
	rm -f binary

-include $(ddMd_tests_chemistry_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

