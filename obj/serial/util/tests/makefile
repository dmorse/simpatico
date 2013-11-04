OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/util/patterns.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/util/tests/sources.mk

TEST=util/tests/Test

all: $(util_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(util_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(OBJ_DIR)/$(TEST) $(SRC_DIR)/util/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in parallel/util/tests/log" > count
	@cat count

clean:
	rm -f $(util_tests_OBJS) $(util_tests_OBJS:.o=.d)
	rm -f $(OBJ_DIR)/$(TEST) 
	rm -f dummy log count

$(OBJ_DIR)/$(TEST): $(OBJ_DIR)/$(TEST).o $(util_LIB)

-include $(util_tests_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)

