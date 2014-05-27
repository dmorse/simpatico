OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/mdCf/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/mdCf/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mdCf/sources.mk
include $(SRC_DIR)/mdCf/tests/sources.mk

all: $(mdCf_tests_OBJS)

run:
	cd chemistry; ./Test
	cd storage; ./Test
	cd configIos; ./Test
	cd processor; ./Test

clean:
	cd chemistry; make clean
	cd storage; make clean
	cd configIos; make clean
	cd processor; make clean

-include $(mdCf_tests_OBJS:.o=.d)
-include $(mdCf_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

