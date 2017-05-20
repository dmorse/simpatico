BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/simp/config.mk
include $(BLD_DIR)/util/config.mk
include $(BLD_DIR)/tools/config.mk
include $(SRC_DIR)/tools/patterns.mk
include $(SRC_DIR)/simp/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/tools/sources.mk
include $(SRC_DIR)/tools/tests/chemistry/sources.mk

all: $(tools_tests_chemistry_OBJS)

clean:
	rm -f $(tools_tests_chemistry_OBJS) 
	rm -f $(tools_tests_chemistry_OBJS:.o=.d)
	rm -f $(tools_tests_chemistry_OBJS:.o=)

-include $(tools_tests_chemistry_OBJS:.o=.d)
-include $(tools_OBJS:.o=.d)
-include $(simp_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

