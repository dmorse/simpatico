OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/ddMd/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/ddMd/patterns.mk
include $(SRC_DIR)/ddMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/ddMd/tests/sources.mk

TEST=ddMd/tests/Test

all: $(ddMd_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(ddMd_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(MPIRUN) 6 $(OBJ_DIR)/$(TEST) $(SRC_DIR)/ddMd/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in src/ddMd/tests/log" > count
	@cat count

clean:
	rm -f $(ddMd_tests_OBJS) $(ddMd_tests_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST)
	rm -f log count
        
-include $(ddMd_tests_OBJS:.o=.d)
-include $(ddMd_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

