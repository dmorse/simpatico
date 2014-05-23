OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/mdPp/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/mdPp/patterns.mk
include $(SRC_DIR)/mdPp/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mdPp/tests/chemistry/sources.mk

all: $(mdPp_tests_chemistry_OBJS)

clean:
	rm -f $(mdPp_tests_chemistry_OBJS) 
	rm -f $(mdPp_tests_chemistry_OBJS:.o=.d)
	rm -f $(mdPp_tests_chemistry_OBJS:.o=)

-include $(mdPp_tests_chemistry_OBJS:.o=.d)
-include $(mdPp_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

