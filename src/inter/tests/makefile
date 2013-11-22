OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/config.mk
include $(OBJ_DIR)/inter/config.mk
include $(OBJ_DIR)/util/config.mk
include $(SRC_DIR)/inter/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/inter/tests/sources.mk

TEST=inter/tests/Test

all: $(inter_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(inter_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(OBJ_DIR)/$(TEST) $(SRC_DIR)/inter/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in inter/tests/log" > count
	@cat count

clean:
	rm -f $(inter_tests_OBJS) $(inter_tests_OBJS:.o=.d) $(OBJ_DIR)/$(TEST)
	rm -f log count
	cd pair; make clean
	cd bond; make clean
	cd angle; make clean
	cd dihedral; make clean

$(OBJ_DIR)/$(TEST): $(OBJ_DIR)/$(TEST).o $(LIBS)

-include $(inter_tests_OBJS:.o=.d)
