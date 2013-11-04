OBJ_DIR_REL =../..
include $(OBJ_DIR_REL)/compiler.mk
include $(OBJ_DIR)/mcMd/defines.mk
include $(OBJ_DIR)/inter/defines.mk
include $(OBJ_DIR)/util/defines.mk
include $(SRC_DIR)/mcMd/patterns.mk
include $(SRC_DIR)/mcMd/sources.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/mcMd/tests/sources.mk

TEST=mcMd/tests/Test

all: $(mcMd_tests_OBJS) $(OBJ_DIR)/$(TEST)

run: $(mcMd_tests_OBJS) $(OBJ_DIR)/$(TEST)
	$(OBJ_DIR)/$(TEST) $(SRC_DIR)/mcMd/tests/ > log
	@echo `grep failed log` ", "\
             `grep successful log` "in src/mcMd/tests/log" > count
	@cat count

clean:
	rm -f $(mcMd_tests_OBJS) $(OBJ_DIR)/$(TEST) $(mcMd_tests_OBJS:.o=.d)
	rm -f *.cfg log count
	cd chemistry; make clean
	cd mcSimulation; make clean
	cd mdSimulation; make clean
	cd neighbor; make clean
	cd species; make clean

$(OBJ_DIR)/$(TEST): $(OBJ_DIR)/$(TEST).o $(LIBS)

-include $(mcMd_tests_OBJS:.o=.d)

