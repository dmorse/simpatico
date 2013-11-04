OBJ_DIR_REL =../../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/mcMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/mcSimulation/sources.mk

TEST=mcMd/tests/mcSimulation/Test

all: $(mcMd_tests_mcSimulation_OBJS) $(OBJ_DIR)/$(TEST)

clean:
	rm -f $(mcMd_tests_mcSimulation_OBJS) $(mcMd_tests_mcSimulation_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 
	rm -f out/param

-include $(mcMd_tests_mcSimulation_OBJS:.o=.d)
-include $(mcMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

